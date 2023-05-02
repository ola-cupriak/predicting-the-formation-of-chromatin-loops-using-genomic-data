"""
This is a boilerplate pipeline 'model_training'
generated using Kedro 0.18.7
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import optimize_parameters
from .nodes import read_data
from .nodes import split_data, save_split_idxes
from .nodes import train_and_eval



def create_pipeline(type: str, **kwargs) -> Pipeline:
    namespace = type
    pipeline_template = pipeline([
        node(
            func=read_data,
            inputs=["concatenated_combined_functional_genomics_data", "params:cell_types", "params:type",
                    "params:features_include_only", "params:features_exclude"],
            outputs="read_data",
            name="read_data_to_dict_node",
        ),
        node(
            func=split_data,
            inputs=["read_data", "params:type", "params:test_size", "params:stratify", "params:random_state"],
            outputs="split_data",
            name="split_data_node",
        ),
        node(
            func=save_split_idxes,
            inputs="split_data",
            outputs="split_data_idxes",
            name="split_data_idx_node",
        ),
        node(
            func=optimize_parameters,
            inputs=["split_data", "params:log_reg.type", "params:log_reg.params", "params:log_reg.optimize", "params:log_reg.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:log_reg.optim_time", "params:log_reg.n_trials"],
            outputs=["logistic_regression_params", "logistic_regression_optimalization_plot"],
            name="optimize_parameters_logistic_regression_node",
        ),
        node(
            func=optimize_parameters,
            inputs=["split_data", "params:rf.type", "params:rf.params", "params:rf.optimize", "params:rf.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:rf.optim_time", "params:rf.n_trials"],
            outputs=["random_forest_params", "random_forest_optimalization_plot"],
            name="optimize_parameters_random_forest_node",
        ),
        node(
            func=optimize_parameters,
            inputs=["split_data", "params:lgbm.type", "params:lgbm.params", "params:lgbm.optimize", "params:lgbm.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:lgbm.optim_time", "params:lgbm.n_trials"],
            outputs=["lightgbm_params", "lightgbm_optimalization_plot"],
            name="optimize_parameters_lightgbm_node",
        ),
        node(
            func=optimize_parameters,
            inputs=["split_data", "params:xgb.type", "params:xgb.params", "params:xgb.optimize", "params:xgb.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:xgb.optim_time", "params:xgb.n_trials"],
            outputs=["xgboost_params", "xgboost_optimalization_plot"],
            name="optimize_parameters_xgboost_node",
        ),
         node(
            func=optimize_parameters,
            inputs=["split_data", "params:dt.type", "params:dt.params", "params:dt.optimize", "params:dt.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:dt.optim_time", "params:dt.n_trials"],
            outputs=["decision_tree_params", "decision_tree_optimalization_plot"],
            name="optimize_parameters_decision_tree_node",
        ),
        node(
            func=train_and_eval,
            inputs=["split_data", "params:log_reg.type", "logistic_regression_params", "params:log_reg.run", "params:run_name"],
            outputs=["logistic_regression_models", "logistic_regression_metrics", "logistic_regression_confusionmatrix",
                     "logistic_regression_feature_importance_df", "logistic_regression_feature_importance_plot"],
            name="train_and_eval_logistic_regression_node",
        ),
        node(
            func=train_and_eval,
            inputs=["split_data", "params:rf.type", "random_forest_params", "params:rf.run", "params:run_name"],
            outputs=["random_forest_models", "random_forest_metrics", "random_forest_confusionmatrix",
                     "random-forest_feature_importance_df", "random-forest_feature_importance_plot"],
            name="train_and_eval_random_forest_node",
        ),
        node(
            func=train_and_eval,
            inputs=["split_data", "params:lgbm.type", "lightgbm_params", "params:lgbm.run", "params:run_name"],
            outputs=["lightgbm_models", "lightgbm_metrics", "lightgbm_confusionmatrix",
                     "lightgbm_feature_importance_df", "lightgbm_feature_importance_plot"],
            name="train_and_eval_lightgbm_node",
        ),
        node(
            func=train_and_eval,
            inputs=["split_data", "params:xgb.type", "xgboost_params", "params:xgb.run", "params:run_name"],
            outputs=["xgboost_models", "xgboost_metrics", "xgboost_confusionmatrix",
                     "xgboost_feature_importance_df", "xgboost_feature_importance_plot"],
            name="train_and_eval_xgboost_node",
        ),
        node(
            func=train_and_eval,
            inputs=["split_data", "params:dt.type", "decision_tree_params", "params:dt.run", "params:run_name"],
            outputs=["decision_tree_models", "decision_tree_metrics", "decision_tree_confusionmatrix",
                     "decision_tree_feature_importance_df", "decision_tree_feature_importance_plot"],
            name="train_and_eval_decision_tree_node",
        ),
    ])

    main_pipeline = pipeline(
        pipe=pipeline_template,
        inputs=["concatenated_combined_functional_genomics_data",],
        parameters=[
            "params:run_name",
            "params:cell_types",
            "params:features_include_only",
            "params:features_exclude",
            "params:test_size",
            "params:validation_size",
            "params:random_state",
            "params:log_reg.type",
            "params:log_reg.params",
            "params:log_reg.run",
            "params:log_reg.optimize",
            "params:log_reg.optim_time",
            "params:log_reg.n_trials",
            "params:rf.type",
            "params:rf.params",
            "params:rf.run",
            "params:rf.optimize",
            "params:rf.optim_time",
            "params:rf.n_trials",
            "params:lgbm.type",
            "params:lgbm.params",
            "params:lgbm.run",
            "params:lgbm.optimize",
            "params:lgbm.optim_time",
            "params:lgbm.n_trials",
            "params:xgb.type",
            "params:xgb.params",
            "params:xgb.run",
            "params:xgb.optimize",
            "params:xgb.optim_time",
            "params:xgb.n_trials",
            "params:dt.type",
            "params:dt.params",
            "params:dt.run",
            "params:dt.optimize",
            "params:dt.optim_time",
            "params:dt.n_trials",
        ],
        namespace=namespace,
    )

    return main_pipeline
