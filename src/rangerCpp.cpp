/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <vector>
#include <sstream>
#include <memory>
#include <utility>

#include "globals.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"
#include "ForestSurvival.h"
#include "ForestProbability.h"
#include "Data.h"
#include "DataChar.h"
#include "DataRcpp.h"
#include "DataFloat.h"
#include "utility.h"

using namespace sirus;

// [[Rcpp::export]]
Rcpp::List rangerCpp(uint treetype, std::string dependent_variable_name, Rcpp::NumericMatrix& input_data,
    std::vector<std::string> variable_names, uint mtry, uint num_trees, bool verbose, uint seed, uint num_threads,
    bool write_forest, uint importance_mode_r, uint min_node_size,
    std::vector<std::vector<double>>& split_select_weights, bool use_split_select_weights,
    std::vector<std::string>& always_split_variable_names, bool use_always_split_variable_names,
    std::string status_variable_name, bool prediction_mode, Rcpp::List loaded_forest, Rcpp::RawMatrix snp_data,
    bool sample_with_replacement, bool probability, std::vector<std::string>& unordered_variable_names,
    bool use_unordered_variable_names, bool save_memory, uint splitrule_r, std::vector<double>& case_weights,
    bool use_case_weights, std::vector<double>& class_weights, bool predict_all, bool keep_inbag,
    std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout, uint prediction_type_r,
    uint num_random_splits, bool order_snps, 
    bool oob_error, uint max_depth, std::vector<std::vector<size_t>>& inbag, bool use_inbag) {

  Rcpp::List result;

  try {
    std::unique_ptr<Forest> forest { };
    std::unique_ptr<Data> data { };

    // Empty split select weights and always split variables if not used
    if (!use_split_select_weights) {
      split_select_weights.clear();
    }
    if (!use_always_split_variable_names) {
      always_split_variable_names.clear();
    }
    if (!use_unordered_variable_names) {
      unordered_variable_names.clear();
    }
    if (!use_case_weights) {
      case_weights.clear();
    }
    if (!use_inbag) {
      inbag.clear();
    }

    std::ostream* verbose_out;
    if (verbose) {
      verbose_out = &Rcpp::Rcout;
    } else {
      verbose_out = new std::stringstream;
    }

    size_t num_rows;
    size_t num_cols;
    num_rows = input_data.nrow();
    num_cols = input_data.ncol();

    // Initialize data 
    data = make_unique<DataRcpp>(input_data, variable_names, num_rows,
      num_cols);

    // If there is snp data, add it
    if (snp_data.nrow() > 1) {
      data->addSnpData(snp_data.begin(), snp_data.ncol());

      // Load SNP order if available
      if (prediction_mode && loaded_forest.containsElementNamed("snp.order")) {
        std::vector<std::vector<size_t>> snp_order = loaded_forest["snp.order"];
        data->setSnpOrder(snp_order);
      }
    }

    switch (treetype) {
    case TREE_CLASSIFICATION:
      if (probability) {
        forest = make_unique<ForestProbability>();
      } else {
        forest = make_unique<ForestClassification>();
      }
      break;
    case TREE_REGRESSION:
      forest = make_unique<ForestRegression>();
      break;
    case TREE_SURVIVAL:
      forest = make_unique<ForestSurvival>();
      break;
    case TREE_PROBABILITY:
      forest = make_unique<ForestProbability>();
      break;
    }

    ImportanceMode importance_mode = (ImportanceMode) importance_mode_r;
    SplitRule splitrule = (SplitRule) splitrule_r;
    PredictionType prediction_type = (PredictionType) prediction_type_r;

    // Init Ranger
    forest->initR(dependent_variable_name, std::move(data), mtry, num_trees, verbose_out, seed, num_threads,
        importance_mode, min_node_size, split_select_weights, always_split_variable_names, status_variable_name,
        prediction_mode, sample_with_replacement, unordered_variable_names, save_memory, splitrule, case_weights,
        inbag, predict_all, keep_inbag, sample_fraction, alpha, minprop, holdout, prediction_type, num_random_splits, 
        order_snps, max_depth);

    // Load forest object if in prediction mode
    if (prediction_mode) {
      size_t dependent_varID = loaded_forest["dependent.varID"];
      //size_t num_trees = loaded_forest["num.trees"];
      std::vector<std::vector<std::vector<size_t>> > child_nodeIDs = loaded_forest["child.nodeIDs"];
      std::vector<std::vector<size_t>> split_varIDs = loaded_forest["split.varIDs"];
      std::vector<std::vector<double>> split_values = loaded_forest["split.values"];
      std::vector<bool> is_ordered = loaded_forest["is.ordered"];

      if (treetype == TREE_CLASSIFICATION) {
        std::vector<double> class_values = loaded_forest["class.values"];
        auto& temp = dynamic_cast<ForestClassification&>(*forest);
        temp.loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values, class_values,
            is_ordered);
      } else if (treetype == TREE_REGRESSION) {
        auto& temp = dynamic_cast<ForestRegression&>(*forest);
        temp.loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values, is_ordered);
      } else if (treetype == TREE_SURVIVAL) {
        size_t status_varID = loaded_forest["status.varID"];
        std::vector<std::vector<std::vector<double>> > chf = loaded_forest["chf"];
        std::vector<double> unique_timepoints = loaded_forest["unique.death.times"];
        auto& temp = dynamic_cast<ForestSurvival&>(*forest);
        temp.loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values, status_varID, chf,
            unique_timepoints, is_ordered);
      } else if (treetype == TREE_PROBABILITY) {
        std::vector<double> class_values = loaded_forest["class.values"];
        std::vector<std::vector<std::vector<double>>> terminal_class_counts = loaded_forest["terminal.class.counts"];
        auto& temp = dynamic_cast<ForestProbability&>(*forest);
        temp.loadForest(dependent_varID, num_trees, child_nodeIDs, split_varIDs, split_values, class_values,
            terminal_class_counts, is_ordered);
      }
    } else {
      // Set class weights
      if (treetype == TREE_CLASSIFICATION && !class_weights.empty()) {
        auto& temp = dynamic_cast<ForestClassification&>(*forest);
        temp.setClassWeights(class_weights);
      } else if (treetype == TREE_PROBABILITY && !class_weights.empty()) {
        auto& temp = dynamic_cast<ForestProbability&>(*forest);
        temp.setClassWeights(class_weights);
      }
    }

    // Run Ranger
    forest->run(false, oob_error);

    if (use_split_select_weights && importance_mode != IMP_NONE) {
      if (verbose_out) {
        *verbose_out
            << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
            << std::endl;
      }
    }

    // Use first non-empty dimension of predictions
    const std::vector<std::vector<std::vector<double>>>& predictions = forest->getPredictions();
    if (predictions.size() == 1) {
      if (predictions[0].size() == 1) {
        result.push_back(forest->getPredictions()[0][0], "predictions");
      } else {
        result.push_back(forest->getPredictions()[0], "predictions");
      }
    } else {
      result.push_back(forest->getPredictions(), "predictions");
    }

    // Return output
    result.push_back(forest->getNumTrees(), "num.trees");
    result.push_back(forest->getNumIndependentVariables(), "num.independent.variables");
    if (treetype == TREE_SURVIVAL) {
      auto& temp = dynamic_cast<ForestSurvival&>(*forest);
      result.push_back(temp.getUniqueTimepoints(), "unique.death.times");
    }
    if (!prediction_mode) {
      result.push_back(forest->getMtry(), "mtry");
      result.push_back(forest->getMinNodeSize(), "min.node.size");
      if (importance_mode != IMP_NONE) {
        result.push_back(forest->getVariableImportance(), "variable.importance");
      }
      result.push_back(forest->getOverallPredictionError(), "prediction.error");
    }
    result.push_back(forest->getForestPaths(), "paths");
    result.push_back(forest->getPathsProba(), "paths.proba");
    if (keep_inbag) {
      result.push_back(forest->getInbagCounts(), "inbag.counts");
    }

    // Save forest if needed
    if (write_forest) {
      Rcpp::List forest_object;
      forest_object.push_back(forest->getDependentVarId(), "dependent.varID");
      forest_object.push_back(forest->getNumTrees(), "num.trees");
      forest_object.push_back(forest->getChildNodeIDs(), "child.nodeIDs");
      forest_object.push_back(forest->getSplitVarIDs(), "split.varIDs");
      forest_object.push_back(forest->getSplitValues(), "split.values");
      forest_object.push_back(forest->getIsOrderedVariable(), "is.ordered");

      if (snp_data.nrow() > 1 && order_snps) {
        // Exclude permuted SNPs (if any)
        std::vector<std::vector<size_t>> snp_order = forest->getSnpOrder();
        forest_object.push_back(std::vector<std::vector<size_t>>(snp_order.begin(), snp_order.begin() + snp_data.ncol()), "snp.order");
      }
      
      if (treetype == TREE_CLASSIFICATION) {
        auto& temp = dynamic_cast<ForestClassification&>(*forest);
        forest_object.push_back(temp.getClassValues(), "class.values");
      } else if (treetype == TREE_PROBABILITY) {
        auto& temp = dynamic_cast<ForestProbability&>(*forest);
        forest_object.push_back(temp.getClassValues(), "class.values");
        forest_object.push_back(temp.getTerminalClassCounts(), "terminal.class.counts");
      } else if (treetype == TREE_SURVIVAL) {
        auto& temp = dynamic_cast<ForestSurvival&>(*forest);
        forest_object.push_back(temp.getStatusVarId(), "status.varID");
        forest_object.push_back(temp.getChf(), "chf");
        forest_object.push_back(temp.getUniqueTimepoints(), "unique.death.times");
      }
      result.push_back(forest_object, "forest");
    }
    
    if (!verbose) {
      delete verbose_out;
    }
  } catch (std::exception& e) {
    if (strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << "Error: " << e.what() << " Ranger will EXIT now." << std::endl;
    }
    return result;
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::List rangerMergeCpp(int numTrees, std::vector<std::vector<std::vector<double> > > paths1, std::vector<int> proba1,
                          std::vector<std::vector<std::vector<double> > > paths2, std::vector<int> proba2){
  Rcpp::List result;

  // merge forests
  std::map<std::vector<std::vector<double> >, int> paths_merge;
  size_t counter1 = 0;
  for (auto& path : paths1) {
    paths_merge[path] = proba1[counter1];
    counter1 += 1;
  }
  std::map<std::vector<std::vector<double> >, int>::iterator it;
  size_t counter2 = 0;
  for (auto& path : paths2) {
    it = paths_merge.find(path);
    if (it == paths_merge.end()){
      paths_merge[path] = proba2[counter2];
    }else{
      paths_merge[path] += proba2[counter2];
    }
    counter2 += 1;
  }
  std::vector<std::vector<std::vector<double>>> result_path;
  std::vector<int> result_proba;
  for (std::map<std::vector<std::vector<double>>, int>::iterator it4 = paths_merge.begin(); it4 != paths_merge.end(); ++it4){
    result_path.push_back(it4->first);
    result_proba.push_back(it4->second);
  }
  result.push_back(result_path);
  result.push_back(result_proba);
  
  // stability metric
  std::vector<double> proba;
  for (auto& prob : result_proba){
    double p = ((double)prob)/numTrees;
    proba.push_back(p);
  }
  std::vector<double> proba_cln = proba;
  std::sort(proba_cln.begin(), proba_cln.end());
  proba_cln.erase(std::unique(proba_cln.begin(), proba_cln.end()), proba_cln.end());
  size_t p0SeqSize = 51;
  if (p0SeqSize > proba_cln.size()){
    p0SeqSize = proba_cln.size();
  }
  std::vector<double> proba100(proba_cln.end() - p0SeqSize, proba_cln.end());
  std::vector<double> p0Seq;
  for (size_t i = 0; i < (proba100.size() - 1); i++) {
    double p0 = (proba100[i] + proba100[i + 1])/2;
    p0Seq.push_back(p0);
  }
  std::vector<double> res;
  double sigma;
  double x;
  for (auto& p0 : p0Seq){
    std::vector<double> phi;
    for (auto& prob : proba) {
      sigma = std::sqrt(prob*(1 - prob)/numTrees);
      x = (p0 - prob)/(sigma*std::sqrt(2));
      x = (1 - std::erf(x))/2;
      phi.push_back(x);
    }
    double numerator = 0;
    double denominator = 0;
    for (auto& x: phi){
      numerator += std::pow(x, 2);
      denominator += x;
    }
    double frac = numerator/denominator;
    res.push_back(frac);
  }
  double result_stability = std::accumulate(res.begin(), res.end(), 0.0)/res.size();
  result.push_back(result_stability);

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector stabilityMetricCpp(int numTrees, std::vector<double> proba){
  
  Rcpp::NumericVector result;
  
  std::vector<double> proba_cln = proba;
  std::sort(proba_cln.begin(), proba_cln.end());
  proba_cln.erase(std::unique(proba_cln.begin(), proba_cln.end()), proba_cln.end());
  size_t p0SeqSize = 51;
  if (p0SeqSize > proba_cln.size()){
    p0SeqSize = proba_cln.size();
  }
  std::vector<double> proba100(proba_cln.end() - p0SeqSize, proba_cln.end());
  std::vector<double> p0Seq;
  for (size_t i = 0; i < (proba100.size() - 1); i++) {
    double p0 = (proba100[i] + proba100[i + 1])/2;
    p0Seq.push_back(p0);
  }
  
  std::vector<double> res;
  double sigma;
  double x;
  
  for (auto& p0 : p0Seq){
    
    std::vector<double> phi;
    for (auto& prob : proba) {
      sigma = std::sqrt(prob*(1 - prob)/numTrees);
      x = (p0 - prob)/(sigma*std::sqrt(2));
      x = (1 - std::erf(x))/2;
      phi.push_back(x);
    }
    
    double numerator = 0;
    double denominator = 0;
    for (auto& x: phi){
      numerator += std::pow(x, 2);
      denominator += x;
    }
    
    double frac = numerator/denominator;
    res.push_back(frac);
    
  }
  
  result = std::accumulate(res.begin(), res.end(), 0.0)/res.size();

  return result;
}
  

