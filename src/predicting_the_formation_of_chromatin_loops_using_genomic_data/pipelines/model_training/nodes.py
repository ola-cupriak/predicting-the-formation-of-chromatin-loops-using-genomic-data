import json
import lightgbm as lgb
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import optuna
import pandas as pd
import plotly.graph_objects as go
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import metrics  
from typing import Any, Dict, List, Tuple, Callable
import xgboost as xgb



def read_data(df: pd.DataFrame, cell_types: list, type: str, 
              features_include_only: list = [], features_exclude: list = []
              ) -> Dict[str, pd.DataFrame]:
    """Reads the data from the data frame and returns the data frame with the
    selected cell types, removes the columns that are not needed for the model.
    Args:
        df: Data frame with concatenated data for model training.
        cell_types: List of cell types to be used for model training.
        type: The type of the model to be trained. Either 'within' or 'across'.
        features_include_only: List of columns to be included in the data frame.
        features_exclude: List of columns to be excluded from the data frame.
    Returns:
        Dictionary with DataFrames with the selected cell types 
        and dropped columns that are not needed for the models.
    """
    if features_include_only and features_exclude:
        print('Warning! Both parameters features_include_only and features_exclude are given. Only the features_include_only parameter will be used.')

    df = df[df['cell_type'].isin(cell_types)]
    df = df.drop(['chr', 'x', 'x_start', 'x_end', 'y', 'y_start', 'y_end'], axis=1)
    if features_include_only:
        if 'label' not in features_include_only:
            features_include_only.append('label')
            features_include_only.append('cell_type')
            assert len(set(df.columns).intersection(set(features_include_only))) == len(features_include_only), print('The features_include_only parameter contains columns that are not in the data frame!')
        df = df.loc[:, features_include_only]
    elif features_exclude:
        assert 'label' not in features_exclude, print('The label column cannot be excluded!')
        assert 'cell_type' not in features_exclude, print('The cell_type column cannot be excluded!')
        assert len(set(df.columns).intersection(set(features_exclude))) == len(features_exclude), print('The features_exclude parameter contains columns that are not in the data frame!')
        df = df.drop(features_exclude, axis=1)
    
    if type == 'within':
        df = df.groupby('cell_type')
        return {cell_type: df_cell for cell_type, df_cell in df}
    elif type == 'across':
        new_df = pd.DataFrame()
        for cell_type in pd.unique(df['cell_type']):
            df_cell_pos = df[(df['cell_type'] == cell_type)&(df['label'] == 1)].sample(frac=0.5, random_state=42, replace=False)
            df_cell_neg = df[(df['cell_type'] == cell_type)&(df['label'] == 0)].sample(frac=0.5, random_state=42, replace=False)
            new_df = pd.concat([new_df, df_cell_pos, df_cell_neg])  
        return {'df': new_df}




def split_data(dfs_dict: Dict[str, pd.DataFrame], 
                type: str, 
                test_size: float = 0.2, 
                stratify: list = [], 
                random_state: int = None
                ) -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    """
    Splits the data into training and test sets.
    Args:
        df: Data frame with the data to be split.
        test_size: The proportion of the data to be used for testing.
        stratify: The column name to be used for stratification.
        random_state: The seed to be used for random number generation.
    Returns:
        A dictionary with tuples with the training and test data frames.
    """
    
    for cell_type, df in dfs_dict.items():
        if type == 'within':
            dfs_dict[cell_type] = train_test_split(df, test_size=test_size, stratify=df.loc[:, stratify], random_state=random_state)
        elif type == 'across':
            pass

    return dfs_dict


def save_split_idxes(dfs_dict: Dict[str, pd.DataFrame],
                     type: str, 
                     ) -> Dict[str, pd.DataFrame]:
    """
    Check indexes in trainset and testset.
    Args:
        dfs_dict: Dictionary with the training and test data frames.
    Returns:
        Dictionary with the DataFrames with indexes of training and test datasets.
    """
    if type == 'across':
        return {}
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


def _train_logistic_regression(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}) -> LogisticRegression:
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


def _train_random_forest(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}) -> RandomForestClassifier():
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


def _train_xgboost(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}) -> xgb.XGBClassifier():
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


def _train_lightgbm(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}) -> lgb.LGBMClassifier():
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


def _train_decision_tree(X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = {}) -> tree.DecisionTreeClassifier():
    """
    Trains a decision tree model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    model = tree.DecisionTreeClassifier(**params)
    model.fit(X_train, y_train)

    return model


def _train_model(df_dict: Dict[str, pd.DataFrame], 
                mtype: str,
                model_type: str = 'log_reg', 
                params: dict = {}) -> Dict[str, Any]:
    """
    Train choosen model on the training data.
    Args:
        df_dict: Dictionary with the training and test data frames.
        mtype: within or across
        model_type: The type of model to be trained.
        random_state: The seed to be used for random number generation.
    Returns:
        A dictionary with the trained models.
    """
    if mtype=='within':
        cell_types = list(df_dict.keys())
    else:
        cell_types = list(pd.unique(df_dict['df']['cell_type']))

    model_dict = {} 
    for cell_type in cell_types:
        if mtype=='within':
            df_train = df_dict[cell_type][0]
        else:
            df_train = df_dict['df'][df_dict['df']['cell_type'] != cell_type]

        print(f'Training {model_type} for {cell_type}...')
        cell_params = params[cell_type]
        if not cell_params:
            cell_params = {}
    
        X_train = df_train.drop(['label', 'cell_type'], axis=1)
        y_train = df_train['label']

        if model_type == 'logistic_regression':
            model = _train_logistic_regression(X_train, y_train, cell_params)
        elif model_type == 'random_forest':
            model = _train_random_forest(X_train, y_train, cell_params)
        elif model_type == 'XGBoost':
            model = _train_xgboost(X_train, y_train, cell_params)
        elif model_type == 'LightGBM':
            model = _train_lightgbm(X_train, y_train, cell_params)
        elif model_type == 'decision_tree':
            model = _train_decision_tree(X_train, y_train, cell_params)
        else:
            raise ValueError(f'Invalid model type: {model_type}')

        model_dict[cell_type] = model

    return model_dict


def _make_prediction(df_test: pd.DataFrame, model) -> np.array:
    """
    Makes predictions on the test data.
    Args:
        df_test: Data frame with the test data.
        model: The trained model.
    Returns:
        A numpy array with the predictions.
    """

    X_test = df_test.drop(['label', 'cell_type'], axis=1)
    y_pred = model.predict(X_test)

    return y_pred


def _evaluate_model(model_dict: dict, df_dict: Dict[str, pd.DataFrame], mtype: str, model_type: str
                    ) -> Tuple[dict, dict]:
    """
    Evaluates the model on the test data.
    Args:
        model: The trained model.
        df_test: Data frame with the test data.
        model_type: The type of model to be trained.
    Returns:
        A tuple with 2 dictionaries - with the evaluation metrics and the confusion matrices plots.
    """
    plt.ioff()
    
    metrics_dict_all = {}
    matrices_dict = {}
    
    for cell_type in model_dict.keys():
        print(f'Evaluating {model_type} for {cell_type}...')
        if mtype=='within':
            df_test = df_dict[cell_type][1]
        else:
            df_test = df_dict['df'][df_dict['df']['cell_type'] == cell_type]

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

        mlflow.log_metrics({f'{cell_type}/{model_type}/auc': metrics.roc_auc_score(y_test, y_pred)})
        metrics_dict = {k: float(str(v)) for k, v in metrics_dict.items()}
    
    matrices_dict = {k+'_confusionmatrix': v for k, v in matrices_dict.items()}

    return metrics_dict_all, matrices_dict


def _get_feature_importances(models_dict: dict, model_type: str) -> Tuple[dict, dict]:
    """
    Get the feature importances for the trained model.
    Args:
        model: The trained model.
        model_type: The type of model.
    Returns:
        A tuple with two dictionaries - fiest with the feature importances DataFrame 
        and second with plot of top20 feature importances.
    """
    feature_importances_plot_dict = {}
    feature_importances_dict = {}

    for cell_type, model in models_dict.items():
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
        elif isinstance(model, RandomForestClassifier) or isinstance(model, tree.DecisionTreeClassifier):
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


def optimize_parameters(df_dict: Dict[str, pd.DataFrame], 
                        model_type: str = 'log_reg', 
                        params: dict = {},
                        optimize: bool = True,
                        run: bool = True,
                        test_size: float = 0.2,
                        validation_size = 0.3,
                        stratify: list = [],
                        random_state: int = None,
                        optim_time: int = None,
                        n_trials: int = 10,
                        eval_metric: Callable = metrics.roc_auc_score,
                        direction: str = "maximize",
                        ) -> Tuple[dict, dict]:
    if not run:
        return {'empty_model': ''}, {'empty_model': go.Figure()}
    
    if not params: 
        params = {}
    params_dict = {k: params for k in df_dict.keys()}
    if not optimize:
        return params_dict, {'empty_model': go.Figure()}

    plot_dict = {}
    validation_size = validation_size/(1-test_size)
    for cell_type, (cell_df, _) in df_dict.items():
        print(f'Optimizing parameters of {model_type} for {cell_type} cell...')
        df_train, df_val = train_test_split(cell_df, test_size=validation_size, stratify=cell_df.loc[:, stratify], random_state=random_state)
        X_train = df_train.drop(['label', 'cell_type'], axis=1)
        y_train = df_train['label']
        X_val = df_val.drop(['label', 'cell_type'], axis=1)
        y_val = df_val['label']
        
        del(df_train)
        del(df_val)

        def objective(trial):
            cell_params = params_dict[cell_type]
            params_to_opt = {}
            for name, val_dict in cell_params.items():
                if val_dict['type'] == 'categorical':
                    params_to_opt[name] = trial.suggest_categorical(name, val_dict['choices'])
                elif val_dict['type'] == 'int':
                    params_to_opt[name] = trial.suggest_int(name, **{k: val_dict[k] for k in set(list(val_dict.keys())) - set(['type'])})
                elif val_dict['type'] == 'float':
                    params_to_opt[name] = trial.suggest_float(name, **{k: val_dict[k] for k in set(list(val_dict.keys())) - set(['type'])})

            
            if model_type == 'logistic_regression':
                model = _train_logistic_regression(X_train, y_train, params_to_opt)
            elif model_type == 'random_forest':
                model = _train_random_forest(X_train, y_train, params_to_opt)
            elif model_type == 'decision_tree':
                model = _train_decision_tree(X_train, y_train, params_to_opt)
            elif model_type == 'LightGBM':
                model = _train_lightgbm(X_train, y_train, params_to_opt)
            elif model_type == 'XGBoost':
                model = _train_xgboost(X_train, y_train, params_to_opt)
            else:
                raise ValueError(f"Wrong model type: {model_type}")
            
            valid_preds = model.predict(X_val)
            valid_score = eval_metric(y_val, valid_preds)

            return valid_score

        if random_state:
            sampler = optuna.samplers.TPESampler(seed=random_state)
        else:
            sampler = optuna.samplers.TPESampler()
        study = optuna.create_study(direction=direction, sampler=sampler)

        if n_trials and optim_time:
            study.optimize(objective, n_trials=n_trials, timeout=optim_time)
        elif optim_time:
            study.optimize(objective, timeout=optim_time)
        elif n_trials:
            study.optimize(objective, n_trials=n_trials)
        else:
            raise ValueError("You have to specify either n_trials or optim_time")

        best_params = study.best_params

        params_dict[cell_type] = best_params
        fig = optuna.visualization.plot_optimization_history(study)
        plot_dict[cell_type] = fig
    return params_dict, plot_dict



def train_and_eval(df_dict: Dict[str, pd.DataFrame],
                    mtype: str,  
                    model_type: str = 'log_reg', 
                    params: dict = {}, 
                    run: bool = True,
                    run_name: str = None,
                    neg_sampling_type: str = None,
                    ) -> Tuple[dict, dict, dict, dict, dict]:
    """
    Trains and evaluates the model.
    Args:
        df_dict: A dictionary with the train and test data frames.
        mtype: within or across
        model_type: The type of model to be trained.
        params: A dictionary with the model parameters.
        run: If True, the model will be trained and evaluated with MLFlow.
        run_name: The name of the run.
    Returns:
        A tuple with dictionaries with the trained models, the evaluation metrics, the confusion matrices,
        the feature importances and the feature importances plot.
    """
    if not run:
        return {'empty_model': ''}, {'empty_model': ''}, {'empty_model': plt.figure()}, {'empty_model': pd.DataFrame()}, {'empty_model': plt.figure()}
    
    if run_name:
        mlflow.set_tag("mlflow.runName", run_name)

    params = _dict_partitions(params)

    model_dict = _train_model(df_dict, mtype, model_type, params)
    metrics_dict_all, matrices_dict = _evaluate_model(model_dict, df_dict, mtype, model_type)
    feature_importances_dict, feature_importances_plot_dict = _get_feature_importances(model_dict, model_type)

    return model_dict, metrics_dict_all, matrices_dict, feature_importances_dict, feature_importances_plot_dict