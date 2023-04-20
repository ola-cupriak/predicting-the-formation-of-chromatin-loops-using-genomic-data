"""
This is a boilerplate pipeline 'model_training'
generated using Kedro 0.18.7
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_data
from .nodes import split_data, save_split_idxes
from .nodes import train_model
from .nodes import evaluate_model



def create_pipeline(type: str, **kwargs) -> Pipeline:
    namespace = type
    pipeline_template = pipeline([
        node(
            func=read_data,
            inputs=["concatenated_combined_functional_genomics_data", "params:cell_types", "params:type"],
            outputs="read_data",
            name="read_data_to_dict_node",
        ),
        node(
            func=split_data,
            inputs=["read_data", "params:test_size", "params:stratify", "params:random_state"],
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
            func=train_model,
            inputs=["split_data", "params:log_reg.type", "params:log_reg.params", "params:log_reg.run"],
            outputs="logistic_regression_models",
            name="train_logistic_regression_node",
        ),
        node(
            func=evaluate_model,
            inputs=["logistic_regression_models", "split_data", "params:log_reg.type"],
            outputs="logistic_regression_metrics",
            name="eval_logistic_regression_node",
        ),
        node(
            func=train_model,
            inputs=["split_data", "params:rf.type", "params:rf.params", "params:rf.run"],
            outputs="random_forest_models",
            name="train_random_forest_node",
        ),
        node(
            func=evaluate_model,
            inputs=["random_forest_models", "split_data", "params:rf.type"],
            outputs="random_forest_metrics",
            name="eval_random_forest_node",
        ),
        node(
            func=train_model,
            inputs=["split_data", "params:lgbm.type", "params:lgbm.params", "params:lgbm.run"],
            outputs="lightgbm_models",
            name="train_lightgbm_node",
        ),
        node(
            func=evaluate_model,
            inputs=["lightgbm_models", "split_data", "params:lgbm.type"],
            outputs="lightgbm_metrics",
            name="eval_lightgbm_node",
        ),
        node(
            func=train_model,
            inputs=["split_data", "params:xgb.type", "params:xgb.params", "params:xgb.run"],
            outputs="xgboost_models",
            name="train_xgboost_node",
        ),
        node(
            func=evaluate_model,
            inputs=["xgboost_models", "split_data", "params:xgb.type"],
            outputs="xgboost_metrics",
            name="eval_xgboost_node",
        ),
    ])

    main_pipeline = pipeline(
        pipe=pipeline_template,
        inputs=["concatenated_combined_functional_genomics_data",],
        parameters=[
            "params:cell_types",
            "params:test_size",
            "params:random_state",
            "params:log_reg.type",
            "params:log_reg.params",
            "params:log_reg.run",
            "params:rf.type",
            "params:rf.params",
            "params:rf.run",
            "params:lgbm.type",
            "params:lgbm.params",
            "params:lgbm.run",
            "params:xgb.type",
            "params:xgb.params",
            "params:xgb.run",
        ],
        namespace=namespace,
    )

    return main_pipeline
