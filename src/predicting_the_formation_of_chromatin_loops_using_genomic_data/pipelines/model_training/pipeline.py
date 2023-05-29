"""
This is a boilerplate pipeline 'model_training'
generated using Kedro 0.18.7
"""
from kedro.pipeline import Pipeline, node, pipeline
from .nodes import optimize_parameters
from .nodes import read_data
from .nodes import split_data, save_split_idxes
from .nodes import train_and_eval
from kedro.config import ConfigLoader
from kedro.framework.project import settings



def create_pipeline(mtype: str, **kwargs) -> Pipeline:
    namespace = mtype

    conf_path = settings.CONF_SOURCE
    conf_loader = ConfigLoader(conf_source=conf_path)
    parameters = conf_loader["parameters"]
    neg_sampling_type = parameters["neg_sampling_type"]

    organism = parameters["organism"]
    if organism == "human":
        input_data = neg_sampling_type+".concatenated_combined_functional_genomics_data"
        to_add = ""
    elif organism == "fly":
        input_data = neg_sampling_type+".FLY_concatenated_combined_functional_genomics_data"
        to_add = "FLY_"
    pipeline_template = pipeline([
        node(
            func=read_data,
            inputs=[input_data, "params:cell_types", "params:type",
                    "params:features_include_only", "params:features_exclude"],
            outputs=to_add+"read_data",
            name="read_data_to_dict_node",
        ),
        node(
            func=split_data,
            inputs=[to_add+"read_data", "params:type", "params:test_size", "params:stratify", "params:random_state"],
            outputs=to_add+"split_data",
            name="split_data_node",
        ),
        node(
            func=save_split_idxes,
            inputs=[to_add+"split_data", "params:type"],
            outputs=to_add+"split_data_idxes",
            name="split_data_idx_node",
        ),
        node(
            func=optimize_parameters,
            inputs=[to_add+"split_data", "params:log_reg.type", "params:log_reg.params", "params:log_reg.optimize", "params:log_reg.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:log_reg.optim_time", "params:log_reg.n_trials"],
            outputs=[to_add+"logistic_regression_params", to_add+"logistic_regression_optimalization_plot"],
            name="optimize_parameters_logistic_regression_node",
        ),
        node(
            func=optimize_parameters,
            inputs=[to_add+"split_data", "params:rf.type", "params:rf.params", "params:rf.optimize", "params:rf.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:rf.optim_time", "params:rf.n_trials"],
            outputs=[to_add+"random_forest_params", to_add+"random_forest_optimalization_plot"],
            name="optimize_parameters_random_forest_node",
        ),
        node(
            func=optimize_parameters,
            inputs=[to_add+"split_data", "params:lgbm.type", "params:lgbm.params", "params:lgbm.optimize", "params:lgbm.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:lgbm.optim_time", "params:lgbm.n_trials"],
            outputs=[to_add+"lightgbm_params", to_add+"lightgbm_optimalization_plot"],
            name="optimize_parameters_lightgbm_node",
        ),
        node(
            func=optimize_parameters,
            inputs=[to_add+"split_data", "params:xgb.type", "params:xgb.params", "params:xgb.optimize", "params:xgb.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:xgb.optim_time", "params:xgb.n_trials"],
            outputs=[to_add+"xgboost_params", to_add+"xgboost_optimalization_plot"],
            name="optimize_parameters_xgboost_node",
        ),
         node(
            func=optimize_parameters,
            inputs=[to_add+"split_data", "params:dt.type", "params:dt.params", "params:dt.optimize", "params:dt.run", "params:test_size", "params:validation_size",
                    "params:stratify", "params:random_state", "params:dt.optim_time", "params:dt.n_trials"],
            outputs=[to_add+"decision_tree_params", to_add+"decision_tree_optimalization_plot"],
            name="optimize_parameters_decision_tree_node",
        ),
        node(
            func=train_and_eval,
            inputs=[to_add+"split_data", "params:type", "params:log_reg.type", "logistic_regression_params", "params:log_reg.run", "params:run_name", "params:neg_sampling_type"],
            outputs=[to_add+"logistic_regression_models", to_add+"logistic_regression_metrics", to_add+"logistic_regression_confusionmatrix",
                     to_add+"logistic_regression_feature_importance_df", to_add+"logistic_regression_feature_importance_plot"],
            name="train_and_eval_logistic_regression_node",
        ),
        node(
            func=train_and_eval,
            inputs=[to_add+"split_data", "params:type", "params:rf.type", "random_forest_params", "params:rf.run", "params:run_name", "params:neg_sampling_type"],
            outputs=[to_add+"random_forest_models", to_add+"random_forest_metrics", to_add+"random_forest_confusionmatrix",
                     to_add+"random-forest_feature_importance_df", to_add+"random-forest_feature_importance_plot"],
            name="train_and_eval_random_forest_node",
        ),
        node(
            func=train_and_eval,
            inputs=[to_add+"split_data", "params:type", "params:lgbm.type", "lightgbm_params", "params:lgbm.run", "params:run_name", "params:neg_sampling_type"],
            outputs=[to_add+"lightgbm_models", to_add+"lightgbm_metrics", to_add+"lightgbm_confusionmatrix",
                     to_add+"lightgbm_feature_importance_df", to_add+"lightgbm_feature_importance_plot"],
            name="train_and_eval_lightgbm_node",
        ),
        node(
            func=train_and_eval,
            inputs=[to_add+"split_data", "params:type", "params:xgb.type", "xgboost_params", "params:xgb.run", "params:run_name", "params:neg_sampling_type"],
            outputs=[to_add+"xgboost_models", to_add+"xgboost_metrics", to_add+"xgboost_confusionmatrix",
                     to_add+"xgboost_feature_importance_df", to_add+"xgboost_feature_importance_plot"],
            name="train_and_eval_xgboost_node",
        ),
        node(
            func=train_and_eval,
            inputs=[to_add+"split_data", "params:type", "params:dt.type", "decision_tree_params", "params:dt.run", "params:run_name", "params:neg_sampling_type"],
            outputs=[to_add+"decision_tree_models", to_add+"decision_tree_metrics", to_add+"decision_tree_confusionmatrix",
                     to_add+"decision_tree_feature_importance_df", to_add+"decision_tree_feature_importance_plot"],
            name="train_and_eval_decision_tree_node",
        ),
    ])

    main_pipeline = pipeline(
        pipe=pipeline_template,
        inputs=[input_data],
        parameters=[
            "params:run_name",
            "params:neg_sampling_type",
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
