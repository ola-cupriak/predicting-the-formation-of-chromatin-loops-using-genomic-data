import json
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import optuna
import pandas as pd
import plotly.graph_objects as go
from sklearn import metrics  
from sklearn.model_selection import LeaveOneOut, KFold
from typing import Any, Dict, List, Tuple, Callable
from tqdm import tqdm
import random
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.model_training.nodes import _filter_features
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.model_training.nodes import _generate_metrics, _generate_confusion_matrix
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.model_training.nodes import _get_feature_importance_single_cell
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.model_training.nodes import _choose_model
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.model_training.nodes import _optimize_parameters



def fly_read_data(df: pd.DataFrame,
              features_include_only: list = None, 
              features_exclude: list = None
              ) -> pd.DataFrame:
    """
    Reads the data from the data frame and returns the data frame with the
    selected cell types, removes the columns that are not needed for the model.
    Args:
        df: Data frame with concatenated data for model training.
        features_include_only: List of columns to be included in the data frame.
        features_exclude: List of columns to be excluded from the data frame.
    Returns:
        Dictionary with DataFrames with the selected cell types 
        and dropped columns that are not needed for the models.
    """
    features_include_only = features_include_only or []
    features_exclude = features_exclude or []

    df = df.drop(['chr', 'x', 'x_start', 'x_end', 'y', 'y_start', 'y_end'], axis=1)
    
    df = _filter_features(df, features_include_only, features_exclude)
    
    return df


def _train_and_pred_with_CV(df: pd.DataFrame, 
                                model_type: str = 'log_reg', 
                                params: dict = None,
                                cv: int = 5) -> Dict[str, Any]:
    """
    Train choosen model with k-fold cross-validation 
    and return the true and predicted labels.
    Args:
        df: pandas DataFrame with the data for model training.
        model_type: The type of model to be trained.
        params: Dictionary with the parameters for the model.
    Returns:
        pandas DataFrame with the true and predicted labels.
    """
    params = params or {}

    true_and_pred = {'true': [], 'pred': []}

    kfoldcv = KFold(n_splits=cv, shuffle=True, random_state=42)
    print(f'Performing {str(cv)}-fold CV for {model_type} model...')
    for train_index, test_index in tqdm(kfoldcv.split(df)):
    
        X_train, X_test = df.iloc[train_index], df.iloc[test_index]
        X_train, X_test = X_train.drop(['label', 'cell_type'], axis=1), X_test.drop(['label', 'cell_type'], axis=1)
        y_train, y_test = df.iloc[train_index]['label'], df.iloc[test_index]['label']

        model = _choose_model(model_type, X_train, y_train, params)
        
        pred = model.predict(X_test)
        true_and_pred['true'] += list(y_test)
        true_and_pred['pred'] += list(pred)

    true_and_pred = pd.DataFrame(true_and_pred)
    return true_and_pred
 

def _evaluate_CV(true_and_pred: pd.DataFrame,
                    model_type: str
                    ) -> Dict[str, Any]:
    """
    Evaluate the CV results and return the metrics.
    Logs the AUC to MLflow.
    Args:
        true_and_pred: pandas DataFrame with the true and predicted labels.
        model_type: The type of model to be trained.
    Returns:
        Dictionary with the metrics and confusion matrix plot.
    """
    plt.ioff()
    y_test = np.array(true_and_pred['true'])
    y_pred = np.array(true_and_pred['pred'])

    mlflow.log_metrics({f'fly_CNS_L3/{model_type}/auc': metrics.roc_auc_score(y_test, y_pred)})
    
    metrics_dict = _generate_metrics(y_test, y_pred, cell_type='fly_CNS_L3', model_type=model_type)
    confusion_matrix = _generate_confusion_matrix(y_test, y_pred, cell_type='fly_CNS_L3', model_type=model_type)
    
    return metrics_dict, confusion_matrix


def fly_optimize_parameters(df: pd.DataFrame, 
                        model_type: str = 'log_reg', 
                        params: dict = None,
                        optimize: bool = True,
                        run: bool = True,
                        cv: int = 5,
                        random_state: int = None,
                        optim_time: int = None,
                        n_trials: int = 10,
                        eval_metric: Callable = metrics.roc_auc_score,
                        direction: str = "maximize",
                        ) -> tuple:
    """
    Optimize the parameters for the model with optuna and cross-validation.
    Args:
        df: pandas DataFrame with the data for model training.
        model_type: The type of model to be trained.
        params: Dictionary with the parameters for the model.
        optimize: If True, the parameters will be optimized.
        run: If True, the model will be trained with the optimized parameters.
        cv: Number of folds for cross-validation.
        random_state: Random state for cross-validation.
        optim_time: Time in seconds for optimization.
        n_trials: Number of trials for optimization.
        eval_metric: Evaluation metric for optimization.
        direction: Direction of optimization.
    Returns:
        A dictionary with the best parameters.
        A optimization history plot.
    """
    random.seed(random_state)
    params = params or {}

    if not run:
        return {}, go.Figure()
    
    if not optimize:
        return params, go.Figure()
    
    print(f'Optimizing parameters of {model_type}...')
    X = df.drop(['label', 'cell_type'], axis=1)
    y = df['label']

    best_params, fig = _optimize_parameters([X, y], 
                                            model_type, params, 
                                            optim_time, n_trials, 
                                            eval_metric, direction, 
                                            random_state, cross_val=cv)

    return best_params, fig


def fly_train_and_eval(df: pd.DataFrame, 
                    model_type: str = 'log_reg', 
                    params: dict = None, 
                    run: bool = True,
                    cv: int = 5,
                    run_name: str = None,
                    neg_sampling_type: str = None,
                    ) -> Tuple[dict, dict, dict, dict, dict]:
    """
    Trains and evaluates the model.
    Args:
        df: pandas DataFrame with the data for model training.
        model_type: The type of model to be trained.
        params: A dictionary with the model parameters.
        run: If True, the model will be trained and evaluated with MLFlow.
        cv: Number of folds for cross-validation.
        run_name: The name of the run.
        neg_sampling_type: The type of negative sampling.
    Returns:
        A tuple with traind model, dictionary with the evaluation metrics, the confusion matrix,
        the feature importances and the feature importances plot.
    """
    params = params or {}

    if not run:
        return {'empty_model': ''}, {'empty_model': plt.figure()}, {'empty_model': pd.DataFrame()}, {'empty_model': plt.figure()}
    
    if run_name:
        mlflow.set_tag("mlflow.runName", run_name)

    print(f'Evaluating {model_type} model...')
    true_and_pred = _train_and_pred_with_CV(df, model_type, params, cv)
    metrics, confusion_matrix = _evaluate_CV(true_and_pred, model_type)

    print(f'Training {model_type} model on full data...')
    X = df.drop(['label', 'cell_type'], axis=1)
    y = df['label']
    model = _choose_model(model_type, X, y, params)

    feature_importances, feature_importances_plot = _get_feature_importance_single_cell(model, cell_type='fly_CNS_L3', model_type=model_type)
    feature_importances = pd.DataFrame(feature_importances, columns=['importances']).reset_index().rename({'index': 'feature'}, axis=1)

    return model, metrics, confusion_matrix, feature_importances, feature_importances_plot