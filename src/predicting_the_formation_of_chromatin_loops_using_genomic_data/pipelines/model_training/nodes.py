import lightgbm as lgb
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import os
import pandas as pd
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn import metrics  
from sklearn import utils
from typing import Any, Callable, Dict
import xgboost as xgb
import yaml



def read_data(df: pd.DataFrame, cell_types: list, type: str) -> pd.DataFrame:
    """Reads the data from the data frame and returns the data frame with the
    selected cell types, removes the columns that are not needed for the model.
    Args:
        df: Data frame with concatenated data for model training.
        cell_types: List of cell types to be used for model training.
        shuffle: If True, the data frame is shuffled.
    Returns:
        Data frame with the selected cell types and dropped columns that are not needed for the model.
    """
    df = df[df['cell_type'].isin(cell_types)]
    df = df.drop(['chr', 'x', 'x_start', 'x_end', 'y', 'y_start', 'y_end'], axis=1)
    
    if type == 'within':
        df = df.groupby('cell_type')
        return {cell_type: df_cell for cell_type, df_cell in df}
    elif type == 'across':
        return {'all': df}


def split_data(dfs_dict: Dict[str, pd.DataFrame], 
                test_size: float = 0.2, 
                stratify: list = [], 
                random_state: int = None) -> tuple:
    """
    Splits the data into training and test sets.
    Args:
        df: Data frame with the data to be split.
        test_size: The proportion of the data to be used for testing.
        stratify: The column name to be used for stratification.
        random_state: The seed to be used for random number generation.
    Returns:
        A tuple with the training and test data frames.
    """
    for cell_type, df in dfs_dict.items():
        dfs_dict[cell_type] = train_test_split(df, test_size=test_size, stratify=df.loc[:, stratify], random_state=random_state)

    return dfs_dict


def save_split_idxes(dfs_dict: Dict[str, pd.DataFrame]) -> dict:
    pass


def _train_logistic_regression(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}):
    """
    Trains a logistic regression model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    model = LogisticRegression(**params)  
    model.fit(X_train, y_train)

    return model


def _train_random_forest(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}):
    """
    Trains a random forest model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    model = RandomForestClassifier(**params)  
    model.fit(X_train, y_train)

    return model


def _train_xgboost(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}):
    """
    Trains a XGBoost model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    model = xgb.XGBClassifier(objective='binary:logistic', **params)
    model.fit(X_train, y_train)

    return model


def _train_lightgbm(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}):
    """
    Trains a LightGBM model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    model = lgb.LGBMClassifier(**params)
    model.fit(X_train, y_train)

    return model


def train_model(df_dict: Dict[str, pd.DataFrame], 
                model_type: str = 'log_reg', 
                params: dict = {}, 
                run: bool = True) -> tuple:
    """
    Train choosen model on the training data.
    Args:
        df_train: Data frame with the training data.
        model_type: The type of model to be trained.
        random_state: The seed to be used for random number generation.
    Returns:
        A tuple with the trained model and the training data.
    """
    if not run:
        return {'empty_model': ''}
    if not params:
        params={}
    
    model_dict = {} 
    for cell_type, (df_train, _) in df_dict.items():
        print(f'Training {model_type} for {cell_type}...')
    
        X_train = df_train.drop(['label', 'cell_type'], axis=1)
        y_train = df_train['label']

        if model_type == 'logistic_regression':
            model = _train_logistic_regression(X_train, y_train, params)
        elif model_type == 'random_forest':
            model = _train_random_forest(X_train, y_train, params)
        elif model_type == 'XGBoost':
            model = _train_xgboost(X_train, y_train, params)
        elif model_type == 'LightGBM':
            model = _train_lightgbm(X_train, y_train, params)
        else:
            raise ValueError(f'Invalid model type: {model_type}')

        model_dict[cell_type] = model

    return model_dict


def evaluate_model(model_dict: dict, df_dict: Dict[str, pd.DataFrame], model_type: str) -> dict:
    """
    Evaluates the model on the test data.
    Args:
        model: The trained model.
        df_test: Data frame with the test data.
        model_type: The type of model to be trained.
    Returns:
        A dictionary with the evaluation metrics.
    """
    model_dict = _dict_partitions(model_dict)
    
    metrics_dict_all = {}
    
    for cell_type in model_dict.keys():
        if cell_type == 'empty_model':
            continue

        print(f'Evaluating {model_type} for {cell_type}...')

        df_test = df_dict[cell_type][1]
        model = model_dict[cell_type]

        X_test = df_test.drop(['label', 'cell_type'], axis=1)
        y_test = df_test['label']

        y_pred = model.predict(X_test)

        metrics_dict = {
            f'accuracy/{model_type}/{cell_type}': metrics.accuracy_score(y_test, y_pred),
            f'precision/{model_type}/{cell_type}': metrics.precision_score(y_test, y_pred),
            f'recall/{model_type}/{cell_type}': metrics.recall_score(y_test, y_pred),
            f'f1/{model_type}/{cell_type}': metrics.f1_score(y_test, y_pred),
            f'auc/{model_type}/{cell_type}': metrics.roc_auc_score(y_test, y_pred),
            f'MCC/{model_type}/{cell_type}': metrics.matthews_corrcoef(y_test, y_pred),
        }

        metrics_dict = {k: round(v, 3) for k, v in metrics_dict.items()}

        mlflow.log_metrics(metrics_dict)

        yaml_string = yaml.dump(metrics_dict)
        metrics_dict_all[cell_type] = yaml_string
    
    if not metrics_dict_all:
        metrics_dict_all['empty_model'] = ''

    return metrics_dict_all


