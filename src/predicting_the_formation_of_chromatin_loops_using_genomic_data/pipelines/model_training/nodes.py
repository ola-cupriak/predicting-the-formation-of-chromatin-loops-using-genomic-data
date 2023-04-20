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



def read_data(df: pd.DataFrame, cell_types: list, type: str, 
              features_include_only: list = [], features_exclude: list = []) -> pd.DataFrame:
    """Reads the data from the data frame and returns the data frame with the
    selected cell types, removes the columns that are not needed for the model.
    Args:
        df: Data frame with concatenated data for model training.
        cell_types: List of cell types to be used for model training.
        type: The type of the model to be trained. Either 'within' or 'across'.
        features_include_only: List of columns to be included in the data frame.
        features_exclude: List of columns to be excluded from the data frame.
    Returns:
        Data frame with the selected cell types and dropped columns that are not needed for the model.
    """
    if features_include_only and features_exclude:
        print('Warning! Both parameters features_include_only and features_exclude are given. Only the features_include_only parameter will be used.')

    df = df[df['cell_type'].isin(cell_types)]
    df = df.drop(['chr', 'x', 'x_start', 'x_end', 'y', 'y_start', 'y_end'], axis=1)
    if features_include_only:
        if 'label' not in features_include_only:
            features_include_only.append('label')
        df = df.loc[:, features_include_only]
    elif features_exclude:
        assert 'label' not in features_exclude, print('The label column cannot be excluded!')
        df = df.drop(features_exclude, axis=1)
    
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
    """
    Check indexes in trainset and testset
    """
    idx_dict = {}
    for cell_type, (train_df, test_df) in dfs_dict.items():
        train_idx = pd.DataFrame(train_df.index, columns=['idx'])
        train_idx['type'] = 'training'
        test_idx = pd.DataFrame(test_df.index, columns=['idx'])
        test_idx['type'] = 'test'

        idx_df = pd.concat([train_df, test_df])

        assert len(set(list(train_idx['idx'])) & set(list(test_idx['idx']))) == 0, print('Training and test datasets have common elements!')

        idx_dict[cell_type] = idx_df
    
    return idx_dict


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


def _train_model(df_dict: Dict[str, pd.DataFrame], 
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


def _make_prediction(df_test: pd.DataFrame, model) -> np.array:


    X_test = df_test.drop(['label', 'cell_type'], axis=1)
    y_pred = model.predict(X_test)

    return y_pred


def _evaluate_model(model_dict: dict, df_dict: Dict[str, pd.DataFrame], model_type: str) -> dict:
    """
    Evaluates the model on the test data.
    Args:
        model: The trained model.
        df_test: Data frame with the test data.
        model_type: The type of model to be trained.
    Returns:
        A dictionary with the evaluation metrics.
    """
    plt.ioff()
    
    metrics_dict_all = {}
    matrices_dict = {}
    
    for cell_type in model_dict.keys():
        if cell_type == 'empty_model':
            continue

        print(f'Evaluating {model_type} for {cell_type}...')

        df_test = df_dict[cell_type][1]
        model = model_dict[cell_type]
        y_test = df_test['label']
        y_pred = _make_prediction(df_test, model)
    
        metrics_dict = {
            f'{cell_type}/{model_type}/accuracy': metrics.accuracy_score(y_test, y_pred),
            f'{cell_type}/{model_type}/precision': metrics.precision_score(y_test, y_pred),
            f'{cell_type}/{model_type}/recall': metrics.recall_score(y_test, y_pred),
            f'{cell_type}/{model_type}/f1': metrics.f1_score(y_test, y_pred),
            f'{cell_type}/{model_type}/auc': metrics.roc_auc_score(y_test, y_pred),
            f'{cell_type}/{model_type}/MCC': metrics.matthews_corrcoef(y_test, y_pred),
        }
        metrics_dict = {k: round(v, 3) for k, v in metrics_dict.items()}

        disp = metrics.ConfusionMatrixDisplay.from_predictions(y_test, y_pred)
        disp.ax_.set_title(f'Confusion matrix for {cell_type} cell for {model_type} model')
        matrices_dict[cell_type] = disp.figure_

        mlflow.log_metrics(metrics_dict)
        metrics_dict_toyaml = {k: float(str(v)) for k, v in metrics_dict.items()}
        yaml_string = yaml.dump(metrics_dict_toyaml)
        metrics_dict_all[cell_type] = yaml_string
    
    matrices_dict = {k+'_confusionmatrix': v for k, v in matrices_dict.items()}

    if not metrics_dict_all:
        metrics_dict_all['empty_model'] = ''
    if not matrices_dict:
        matrices_dict['empty_model'] = plt.figure()
    

    return metrics_dict_all, matrices_dict


def _get_feature_importances(models_dict: dict, model_type: str) -> dict:
    feature_importances_plot_dict = {}
    feature_importances_dict = {}

    for cell_type, model in models_dict.items():
        if cell_type == 'empty_model':
            continue
        if isinstance(model, xgb.XGBClassifier):
            importances = model.feature_importances_
            idxs = np.argsort(importances)
            importances = importances[idxs]
            names = np.array(model.get_booster().feature_names)[idxs]
        elif isinstance(model, lgb.LGBMClassifier):
            importances = model.feature_importances_
            idxs = np.argsort(importances)
            importances = importances[idxs]
            names = np.array(model.booster_.feature_name())[idxs]
        elif isinstance(model, LogisticRegression):
            importances = model.coef_[0]
            idxs = np.argsort(importances)
            importances = importances[idxs]
            names = model.feature_names_in_[idxs]
        elif isinstance(model, RandomForestClassifier):
            importances = model.feature_importances_
            idxs = np.argsort(importances)
            importances = importances[idxs]
            names = model.feature_names_in_[idxs]

        importances = pd.Series(importances, index=names)
        fig, ax = plt.subplots()
        importances[-20:].plot.bar(ax=ax)
        ax.set_title(f"Top20 feature importances for {cell_type} cell for {model_type} model")
        ax.set_ylabel("Feature importance")
        fig.tight_layout()

        feature_importances_plot_dict[cell_type] = fig
        feature_importances_dict[cell_type] = pd.DataFrame(importances, columns=['importance'])

    return feature_importances_dict, feature_importances_plot_dict



def train_and_eval(df_dict: Dict[str, pd.DataFrame], 
                    model_type: str = 'log_reg', 
                    params: dict = {}, 
                    run: bool = True):
    
    model_dict = _train_model(df_dict, model_type, params, run)
    metrics_dict_all, matrices_dict = _evaluate_model(model_dict, df_dict, model_type)
    feature_importances_dict, feature_importances_plot_dict = _get_feature_importances(model_dict, model_type)

    return model_dict, metrics_dict_all, matrices_dict, feature_importances_dict, feature_importances_plot_dict