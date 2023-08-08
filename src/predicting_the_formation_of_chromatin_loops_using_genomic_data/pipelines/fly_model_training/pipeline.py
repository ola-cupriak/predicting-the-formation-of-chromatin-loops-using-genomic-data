"""
This is a boilerplate pipeline 'model_training'
generated using Kedro 0.18.7
"""
from kedro.pipeline import Pipeline, node, pipeline
from .nodes import fly_optimize_parameters
from .nodes import fly_read_data
from .nodes import fly_train_and_eval
from kedro.config import ConfigLoader
from kedro.framework.project import settings



def create_pipeline(**kwargs) -> Pipeline:

    conf_path = settings.CONF_SOURCE
    conf_loader = ConfigLoader(conf_source=conf_path)
    parameters = conf_loader["parameters"]
    fly_neg_sampling_type = parameters["fly_neg_sampling_type"]

    input_data = fly_neg_sampling_type+".FLY_concatenated_combined_functional_genomics_data"

    main_pipeline = pipeline([
        node(
            func=fly_read_data,
            inputs=[input_data, "params:fly_features_include_only", "params:fly_features_exclude"],
            outputs="fly_read_data",
            name="fly_read_data_to_dict_node",
        ),
        node(
            func=fly_optimize_parameters,
            inputs=["fly_read_data", "params:fly_log_reg.type", "params:fly_log_reg.params", "params:fly_log_reg.optimize", 
                    "params:fly_log_reg.run", "params:fly_log_reg.cv", "params:random_state", "params:fly_log_reg.optim_time", "params:fly_log_reg.n_trials"],
            outputs=["FLY_logistic_regression_params", "FLY_logistic_regression_optimalization_plot"],
            name="fly_optimize_parameters_logistic_regression_node",
        ),
        node(
            func=fly_optimize_parameters,
            inputs=["fly_read_data", "params:fly_rf.type", "params:fly_rf.params", "params:fly_rf.optimize", 
                    "params:fly_rf.run", "params:fly_rf.cv", "params:random_state", "params:fly_rf.optim_time", "params:fly_rf.n_trials"],
            outputs=["FLY_random_forest_params", "FLY_random_forest_optimalization_plot"],
            name="fly_optimize_parameters_random_forest_node",
        ),
        node(
            func=fly_optimize_parameters,
            inputs=["fly_read_data", "params:fly_lgbm.type", "params:fly_lgbm.params", "params:fly_lgbm.optimize", 
                    "params:fly_lgbm.run", "params:fly_lgbm.cv", "params:random_state", "params:fly_lgbm.optim_time", "params:fly_lgbm.n_trials"],
            outputs=["FLY_lightgbm_params", "FLY_lightgbm_optimalization_plot"],
            name="fly_optimize_parameters_lightgbm_node",
        ),
         node(
            func=fly_optimize_parameters,
            inputs=["fly_read_data", "params:fly_dt.type", "params:fly_dt.params", "params:fly_dt.optimize", 
                    "params:fly_dt.run", "params:fly_dt.cv", "params:random_state", "params:fly_dt.optim_time", "params:fly_dt.n_trials"],
            outputs=["FLY_decision_tree_params", "FLY_decision_tree_optimalization_plot"],
            name="fly_optimize_parameters_decision_tree_node",
        ),
        node(
            func=fly_train_and_eval,
            inputs=["fly_read_data", "params:fly_log_reg.type", "FLY_logistic_regression_params", "params:fly_log_reg.run", "params:fly_log_reg.cv", 
                     "params:fly_run_name", "params:fly_neg_sampling_type"],
            outputs=["FLY_logistic_regression_models", "FLY_logistic_regression_metrics", "FLY_logistic_regression_confusionmatrix",
                     "FLY_logistic_regression_roccurve",
                     "FLY_logistic_regression_feature_importance_df", "FLY_logistic_regression_feature_importance_plot"],
            name="fly_train_and_eval_logistic_regression_node",
        ),
        node(
            func=fly_train_and_eval,
            inputs=["fly_read_data", "params:fly_rf.type", "FLY_random_forest_params", "params:fly_rf.run",  "params:fly_rf.cv",
                      "params:fly_run_name", "params:fly_neg_sampling_type"],
            outputs=["FLY_random_forest_models", "FLY_random_forest_metrics", "FLY_random_forest_confusionmatrix",
                        "FLY_random_forest_roccurve",
                     "FLY_random_forest_feature_importance_df", "FLY_random_forest_feature_importance_plot"],
            name="fly_train_and_eval_random_forest_node",
        ),
        node(
            func=fly_train_and_eval,
            inputs=["fly_read_data", "params:fly_lgbm.type", "FLY_lightgbm_params", "params:fly_lgbm.run", "params:fly_lgbm.cv",
                    "params:fly_run_name", "params:fly_neg_sampling_type"],
            outputs=["FLY_lightgbm_models", "FLY_lightgbm_metrics", "FLY_lightgbm_confusionmatrix",
                        "FLY_lightgbm_roccurve",
                     "FLY_lightgbm_feature_importance_df", "FLY_lightgbm_feature_importance_plot"],
            name="fly_train_and_eval_lightgbm_node",
        ),
        node(
            func=fly_train_and_eval,
            inputs=["fly_read_data", "params:fly_dt.type", "FLY_decision_tree_params", "params:fly_dt.run", "params:fly_dt.cv",
                      "params:fly_run_name", "params:fly_neg_sampling_type"],
            outputs=["FLY_decision_tree_models", "FLY_decision_tree_metrics", "FLY_decision_tree_confusionmatrix",
                        "FLY_decision_tree_roccurve",
                     "FLY_decision_tree_feature_importance_df", "FLY_decision_tree_feature_importance_plot"],
            name="fly_train_and_eval_decision_tree_node",
        ),
    ])

    return main_pipeline
