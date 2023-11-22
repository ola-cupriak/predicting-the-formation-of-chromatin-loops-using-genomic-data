import json
import yaml
import os
import lightgbm as lgb
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import optuna
import pandas as pd
import plotly.graph_objects as go
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import (
    _dict_partitions,
)
from sklearn import tree
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from typing import Any, Dict, List, Tuple, Callable
import random
import time
import shap


def _filter_features(
    df: pd.DataFrame, features_include_only: list, features_exclude: list
) -> pd.DataFrame:
    """Filters the features in the data frame.
    Args:
        df: Data frame with the data to be filtered.
        features_include_only: List of columns to be included in the data frame.
        features_exclude: List of columns to be excluded from the data frame.
    Returns:
        Data frame with the filtered features.
    """
    if features_include_only and features_exclude:
        print(
            "Warning! Both parameters features_include_only and features_exclude are given. Only the features_include_only parameter will be used."
        )
    if features_include_only:
        if "label" not in features_include_only:
            features_include_only.append("label")
            features_include_only.append("cell_type")
            assert len(set(df.columns).intersection(set(features_include_only))) == len(
                features_include_only
            ), print(
                "The features_include_only parameter contains columns that are not in the data frame!"
            )
        df = df.loc[:, features_include_only]
    elif features_exclude:
        assert "label" not in features_exclude, print(
            "The label column cannot be excluded!"
        )
        assert "cell_type" not in features_exclude, print(
            "The cell_type column cannot be excluded!"
        )
        assert len(set(df.columns).intersection(set(features_exclude))) == len(
            features_exclude
        ), print(
            "The features_exclude parameter contains columns that are not in the data frame!"
        )
        df = df.drop(features_exclude, axis=1)

    return df


def read_data(
    df: pd.DataFrame,
    cell_types: list,
    mtype: str,
    random_state: int = None,
    features_include_only: list = None,
    features_exclude: list = None,
    data_fraction: float = 1,
) -> Dict[str, pd.DataFrame]:
    """Reads the data from the data frame and returns the data frame with the
    selected cell types, removes the columns that are not needed for the model.
    Args:
        df: Data frame with concatenated data for model training.
        cell_types: List of cell types to be used for model training.
        mtype: The type of the model to be trained. Either 'within' or 'across'.
        random_state: The seed to be used for random number generation.
        features_include_only: List of columns to be included in the data frame.
        features_exclude: List of columns to be excluded from the data frame.
        data_fraction: The fraction of the data to be used for model training.
    Returns:
        Dictionary with DataFrames with the selected cell types 
        and dropped columns that are not needed for the models.
    """
    features_include_only = features_include_only or []
    features_exclude = features_exclude or []

    df = df[df["cell_type"].isin(cell_types)]
    df = df.drop(["chr", "x", "x_start", "x_end", "y", "y_start", "y_end"], axis=1)

    df = _filter_features(df, features_include_only, features_exclude)

    if mtype == "within":
        df = df.groupby("cell_type")
        return {cell_type: df_cell for cell_type, df_cell in df}
    elif mtype == "across":
        if data_fraction < 1:
            new_df = pd.DataFrame()
            for cell_type in pd.unique(df["cell_type"]):
                df_cell_pos = df[
                    (df["cell_type"] == cell_type) & (df["label"] == 1)
                ].sample(frac=data_fraction, random_state=random_state, replace=False)
                df_cell_neg = df[
                    (df["cell_type"] == cell_type) & (df["label"] == 0)
                ].sample(frac=data_fraction, random_state=random_state, replace=False)
                new_df = pd.concat([new_df, df_cell_pos, df_cell_neg])
        else:
            new_df = df    

        return {"df": new_df}


def split_data(
    dfs_dict: Dict[str, pd.DataFrame],
    type: str,
    test_size: float = 0.2,
    stratify: list = None,
    random_state: int = None,
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
    stratify = stratify or []
    for cell_type, df in dfs_dict.items():
        if type == "within":
            dfs_dict[cell_type] = train_test_split(
                df,
                test_size=test_size,
                stratify=df.loc[:, stratify],
                random_state=random_state,
            )
        elif type == "across":
            pass

    return dfs_dict


def save_split_idxes(
    dfs_dict: Dict[str, pd.DataFrame],
    type: str,
) -> Dict[str, pd.DataFrame]:
    """
    Check indexes in trainset and testset.
    Args:
        dfs_dict: Dictionary with the training and test data frames.
    Returns:
        Dictionary with the DataFrames with indexes of training and test datasets.
    """
    if type == "across":
        return {}
    idx_dict = {}
    for cell_type, (train_df, test_df) in dfs_dict.items():
        train_idx = pd.DataFrame(train_df.index, columns=["idx"])
        train_idx["type"] = "training"
        test_idx = pd.DataFrame(test_df.index, columns=["idx"])
        test_idx["type"] = "test"

        idx_df = pd.concat([train_df, test_df])

        assert (
            len(set(list(train_idx["idx"])) & set(list(test_idx["idx"]))) == 0
        ), print("Training and test datasets have common elements!")

        idx_dict[cell_type] = idx_df

    return idx_dict


def _train_logistic_regression(
    X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = None
) -> LogisticRegression:
    """
    Trains a logistic regression model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    params = params or {}
    model = LogisticRegression(**params)
    model.fit(X_train, y_train)

    return model


def _train_random_forest(
    X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = None
) -> RandomForestClassifier:
    """
    Trains a random forest model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    params = params or {}
    model = RandomForestClassifier(**params)
    model.fit(X_train, y_train)

    return model


def _train_lightgbm(
    X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = None
) -> lgb.LGBMClassifier:
    """
    Trains a LightGBM model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    params = params or {}
    model = lgb.LGBMClassifier(**params)
    model.fit(X_train, y_train)

    return model


def _train_decision_tree(
    X_train: pd.DataFrame, y_train: pd.DataFrame, params: dict = None
) -> tree.DecisionTreeClassifier:
    """
    Trains a decision tree model on the training data and returns the trained model.
    Args:
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    params = params or {}
    model = tree.DecisionTreeClassifier(**params)
    model.fit(X_train, y_train)

    return model


def _choose_model(model_type: str, X_train, y_train, params: dict = None) -> Callable:
    """
    Trains a model on the training data and returns the trained model.
    Args:
        model_type: Type of model to be trained.
        X_train: Data frame with the training data.
        y_train: Data frame with the training labels.
        params: Dictionary with the parameters for the model.
    Returns:
        The trained model.
    """
    params = params or {}

    if model_type == "logistic_regression":
        model = _train_logistic_regression(X_train, y_train, params)
    elif model_type == "random_forest":
        model = _train_random_forest(X_train, y_train, params)
    elif model_type == "LightGBM":
        model = _train_lightgbm(X_train, y_train, params)
    elif model_type == "decision_tree":
        model = _train_decision_tree(X_train, y_train, params)
    else:
        raise ValueError(f"Invalid model type: {model_type}")

    return model


def _create_SHAP_plot(model, X_test, cell_type, model_type):
    fig, ax = plt.subplots()
    explainer = shap.Explainer(model)
    shap_values = explainer.shap_values(X_test)
    shap.summary_plot(shap_values[1], X_test, show=False)
    ax.tick_params(axis="y", labelsize=10)
    odm = "trenowanego"
    mtype = "pomiędzy"
    negtype = "B"
    if model_type == "random_forest":
        model_type = "lasu losowego"
    elif model_type == "logistic_regression":
        model_type = "regresji logistycznej"
        odm = "trenowanej"
    elif model_type == "decision_tree":
        model_type = "drzewa decyzyjnego"
    ax.set_title(
        f'Wykres wartości SHAP dla {model_type} typu "{mtype}"\n{odm} na danych wygenerowanych metodą {negtype}',
        fontsize=12,
    )
    ax.set_xlabel("Wartość SHAP (wpływ na predykcję dla klasy pozytywnej)", fontsize=10)
    fig.tight_layout()
    mlflow.log_figure(fig, f"{cell_type}/{model_type}/SHAP_plot.png")


def _train_model(
    df_dict: Dict[str, pd.DataFrame],
    mtype: str,
    model_type: str = "log_reg",
    params: dict = None,
    SHAP: bool = False,
) -> Dict[str, Any]:
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
    params = params or {}

    if mtype == "within":
        cell_types = list(df_dict.keys())
    else:
        cell_types = list(pd.unique(df_dict["df"]["cell_type"]))

    model_dict = {}
    for cell_type in cell_types:
        if mtype == "within":
            df_train = df_dict[cell_type][0]
        else:
            df_train = df_dict["df"][df_dict["df"]["cell_type"] != cell_type]

        print(f"Training {model_type} for {cell_type}...")

        if cell_type in params.keys():
            cell_params = params[cell_type]
        else:
            cell_params = {}

        X_train = df_train.drop(["label", "cell_type"], axis=1)
        y_train = df_train["label"]

        model = _choose_model(model_type, X_train, y_train, cell_params)

        model_dict[cell_type] = model

        shap_cell = input(f"Do you want to generate SHAP plot for cell {cell_type}? (y/n): ")
        if SHAP and shap_cell == "y":
            _create_SHAP_plot(model, X_train, cell_type, model_type)

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
    X_test = df_test.drop(["label", "cell_type"], axis=1)
    y_pred = model.predict(X_test)

    return y_pred


def _generate_metrics(
    y_test: pd.DataFrame, y_pred: np.array, cell_type: str, model_type: str
) -> dict:
    """
    Generates metrics for the predictions.
    Args:
        y_test: Data frame with the test labels.
        y_pred: Numpy array with the predictions.
        cell_type: The cell type for which the metrics are generated.
        model_type: The type of model for which the metrics are generated.
    Returns:
        A dictionary with the metrics (accuracy, precision, recall, f1, auc, MCC).
    """
    metrics_dict = {
        f"{cell_type}/{model_type}/accuracy": metrics.accuracy_score(y_test, y_pred),
        f"{cell_type}/{model_type}/precision": metrics.precision_score(y_test, y_pred),
        f"{cell_type}/{model_type}/recall": metrics.recall_score(y_test, y_pred),
        f"{cell_type}/{model_type}/f1": metrics.f1_score(y_test, y_pred),
        f"{cell_type}/{model_type}/auc": metrics.roc_auc_score(y_test, y_pred),
        f"{cell_type}/{model_type}/MCC": metrics.matthews_corrcoef(y_test, y_pred),
    }
    metrics_dict = {k: round(v, 3) for k, v in metrics_dict.items()}
    metrics_dict = {k: float(str(v)) for k, v in metrics_dict.items()}

    return metrics_dict


def _generate_confusion_matrix(
    y_test: pd.DataFrame, y_pred: np.array, cell_type: str, model_type: str
):
    """
    Generates confusion matrix for the predictions.
    Args:
        y_test: Data frame with the test labels.
        y_pred: Numpy array with the predictions.
        cell_type: The cell type for which the metrics are generated.
        model_type: The type of model for which the metrics are generated.
    Returns:
        A confusion matrix plot.
    """
    disp = metrics.ConfusionMatrixDisplay.from_predictions(y_test, y_pred)
    disp.ax_.set_title(f"Confusion matrix for {cell_type} cell for {model_type} model")

    return disp.figure_


def _evaluate_model(
    model_dict: dict, df_dict: Dict[str, pd.DataFrame], mtype: str, model_type: str
) -> Tuple[dict, dict]:
    """
    Evaluates the model on the test data.
    Log the AUC to MLflow.
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
        print(f"Evaluating {model_type} for {cell_type}...")
        if mtype == "within":
            df_test = df_dict[cell_type][1]
        else:
            df_test = df_dict["df"][df_dict["df"]["cell_type"] == cell_type]

        model = model_dict[cell_type]
        y_test = df_test["label"]
        y_pred = _make_prediction(df_test, model)

        mlflow.log_metrics(
            {f"{cell_type}/{model_type}/auc": metrics.roc_auc_score(y_test, y_pred)}
        )
        mlflow.log_metrics(
            {f"{cell_type}/{model_type}/acc": metrics.accuracy_score(y_test, y_pred)}
        )

        matrices_dict[cell_type] = _generate_confusion_matrix(
            y_test, y_pred, cell_type, model_type
        )
        metrics_dict_all[cell_type] = _generate_metrics(
            y_test, y_pred, cell_type, model_type
        )

    # log mean auc for all cells
    mlflow.log_metrics(
        {
            f"all_cells/{model_type}/auc": np.mean(
                [v[f"{k}/{model_type}/auc"] for k, v in metrics_dict_all.items()]
            )
        }
    )
    mlflow.log_metrics(
        {
            f"all_cells/{model_type}/acc": np.mean(
                [v[f"{k}/{model_type}/accuracy"] for k, v in metrics_dict_all.items()]
            )
        }
    )

    matrices_dict = {k + "_confusionmatrix": v for k, v in matrices_dict.items()}

    return metrics_dict_all, matrices_dict


def _get_feature_importance_single_cell(
    model, cell_type: str, model_type: str
) -> Tuple[dict, dict]:
    """
    Get the feature importances for the trained model for a single cell type.
    Args:
        model: The trained model.
        model_type: The type of model.
    Returns:
        A pandas data frame with the feature importances and a plot with the feature importances.
    """
    if isinstance(model, lgb.LGBMClassifier):
        importances = model.feature_importances_
        idxs = np.argsort(importances)
        importances = importances[idxs]
        names = np.array(model.booster_.feature_name())[idxs]
    elif isinstance(model, LogisticRegression):
        importances = model.coef_[0]
        idxs = np.argsort(importances)
        importances = importances[idxs]
        names = model.feature_names_in_[idxs]
    elif isinstance(model, RandomForestClassifier) or isinstance(
        model, tree.DecisionTreeClassifier
    ):
        importances = model.feature_importances_
        idxs = np.argsort(importances)
        importances = importances[idxs]
        names = model.feature_names_in_[idxs]

    importances = pd.Series(importances, index=names)
    fig, ax = plt.subplots()
    importances[-20:].plot.barh(ax=ax)
    ax.set_title(
        f"Top 20 feature importances for {cell_type} cell for {model_type} model"
    )
    ax.set_xlabel("Feature importance")
    fig.tight_layout()

    return importances, fig


def _get_feature_importances(models_dict: dict, model_type: str) -> Tuple[dict, dict]:
    """
    Get the feature importances for the trained models for all cell types.
    Args:
        models_dict: The dictionary with the trained models for all cell types.
        model_type: The type of model.
    Returns:
        A tuple with two dictionaries - fiest with the feature importances DataFrame
        and second with plot of top20 feature importances.
    """
    feature_importances_plot_dict = {}
    feature_importances_dict = {}

    for cell_type, model in models_dict.items():
        importances, fig = _get_feature_importance_single_cell(
            model, cell_type, model_type
        )
        feature_importances_plot_dict[cell_type] = fig
        feature_importances_dict[cell_type] = pd.DataFrame(
            importances, columns=["importance"]
        )

    return feature_importances_dict, feature_importances_plot_dict


def _train_predic_eval(
    model_type, X_train, y_train, X_val, y_val, eval_metric, params: dict = None
) -> float:
    """
    Trains the model, makes predictions and evaluates the model.
    Args:
        model_type: The type of model to be trained.
        X_train: The features for the training set.
        y_train: The target for the training set.
        X_val: The features for the validation set.
        y_val: The target for the validation set.
        eval_metric: The evaluation metric.
        params: The parameters for the model.
    Returns:
        The validation score.
    """
    params = params or {}
    model = _choose_model(model_type, X_train, y_train, params)
    valid_preds = model.predict(X_val)
    valid_score = eval_metric(y_val, valid_preds)

    return valid_score


def _cross_val(model_type, X, y, folds: int, params: dict = None) -> float:
    """
    Performs cross validation for the given model type.
    Args:
        model_type: The type of model to be trained.
        X: The features.
        y: The target.
        eval_metric: The evaluation metric.
        folds: The number of folds for the cross validation.
        params: The parameters for the model.
    Returns:
        The mean cross validation score.
    """
    params = params or {}
    if model_type == "logistic_regression":
        model = LogisticRegression(**params)
    elif model_type == "random_forest":
        model = RandomForestClassifier(**params)
    elif model_type == "LightGBM":
        model = lgb.LGBMClassifier(**params)
    elif model_type == "decision_tree":
        model = tree.DecisionTreeClassifier(**params)
    else:
        raise ValueError(f"Invalid model type: {model_type}")

    cv_scores = cross_val_score(model, X, y, cv=folds, scoring="roc_auc")

    return cv_scores.mean()


def _optimize_parameters(
    data: list,
    model_type: str,
    params: dict,
    optim_time,
    n_trials,
    eval_metric,
    direction,
    random_state,
    cross_val: int = 0,
) -> tuple:
    """
    Optimize the hyperparameters for a single model.
    Args:
        data: list of datasets for training and validation.
        model_type: The type of model.
        params: The parameters to be optimized.
        optim_time: The time for optimization.
        n_trials: The number of trials.
        eval_metric: The metric to be used for evaluation.
        direction: The direction of optimization.
        random_state: The random state.
        cross_val: The number of folds for cross validation (default: 0 - no cross validation).
    Returns:
        A dictionary with the best parameters.
        A plot with the optimization history.
    """

    def objective(trial):
        params_to_opt = {}
        for name, val_dict in params.items():
            if val_dict["type"] == "categorical":
                params_to_opt[name] = trial.suggest_categorical(
                    name, val_dict["choices"]
                )
            elif val_dict["type"] == "int":
                params_to_opt[name] = trial.suggest_int(
                    name,
                    **{
                        k: val_dict[k]
                        for k in set(list(val_dict.keys())) - set(["type"])
                    },
                )
            elif val_dict["type"] == "float":
                params_to_opt[name] = trial.suggest_float(
                    name,
                    **{
                        k: val_dict[k]
                        for k in set(list(val_dict.keys())) - set(["type"])
                    },
                )
            elif val_dict["type"] == "uniform":
                params_to_opt[name] = trial.suggest_float(
                    name,
                    **{
                        k: val_dict[k]
                        for k in set(list(val_dict.keys())) - set(["type"])
                    },
                )

        if cross_val:
            X, y = data
            valid_score = _cross_val(
                model_type, X, y, folds=cross_val, params=params_to_opt
            )
        else:
            X_train, y_train, X_val, y_val = data
            valid_score = _train_predic_eval(
                model_type, X_train, y_train, X_val, y_val, eval_metric, params_to_opt
            )

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

    # check if the best parameters are better than the default ones
    if cross_val:
        X, y = data
        valid_score_default = _cross_val(model_type, X, y, folds=cross_val)
    else:
        X_train, y_train, X_val, y_val = data
        valid_score_default = _train_predic_eval(
            model_type, X_train, y_train, X_val, y_val, eval_metric
        )

    if valid_score_default > study.best_value:
        print(
            f"Default parameters are better than the optimized ones. Default score: {valid_score_default}, optimized score: {study.best_value}. Using default parameters."
        )
        best_params = {}

    fig = optuna.visualization.plot_optimization_history(study)

    return best_params, fig


def optimize_parameters(
    df_dict: Dict[str, pd.DataFrame],
    mtype: str,
    model_type: str = "log_reg",
    params: dict = None,
    optimize: bool = True,
    run: bool = True,
    test_size: float = 0.2,
    validation_size=0.3,
    stratify: list = None,
    random_state: int = None,
    optim_time: int = None,
    n_trials: int = 10,
    eval_metric: Callable = metrics.roc_auc_score,
    direction: str = "maximize",
) -> tuple:
    """
    Optimize the hyperparameters for the models for all cell types.
    Args:
        df_dict: The dictionary with the data for all cell types.
        mtype: The type of model (within or across).
        model_type: The type of model.
        params: The parameters to be optimized.
        optimize: Whether to optimize the parameters.
        run: Whether to run the optimization.
        test_size: The size of the test set.
        validation_size: The size of the validation set.
        stratify: The column to be used for stratification.
        random_state: The random state.
        optim_time: The time for optimization.
        n_trials: The number of trials.
        eval_metric: The metric to be used for evaluation.
        direction: The direction of optimization.
    Returns:
        A dictionary with the best parameters.
        A optimization history plot.
    """
    random.seed(random_state)
    stratify = stratify or []
    params = params or {}

    if not run or not optimize:
        for mt in ["logistic_regression", "random_forest", "LightGBM", "decision_tree"]:
            if model_type == mt:
                if mt == "LightGBM":
                    mt = "lightgbm"
                if f"{mtype}_{mt}_params.yml" in os.listdir(
                    f"data/05_model_input/Homo_sapiens/{mtype}"
                ):
                    with open(
                        f"data/05_model_input/Homo_sapiens/{mtype}/{mtype}_{mt}_params.yml",
                        "r",
                    ) as f:
                        params = yaml.safe_load(f)
        params = params or {}
        return params, {}

    print(f"Optimizing parameters of {model_type}...")
    params_opt_dict = {}
    optim_fig_dict = {}

    if mtype == "within":
        cells = df_dict.keys()
    else:
        cells = list(set(df_dict["df"]["cell_type"]))

    for cell in cells:
        print(f"Optimizing parameters for cell {cell}...")
        # Choose cell type to validate on
        if mtype == "within":
            validation_size_within = validation_size / (1 - test_size)
            df_train, df_val = train_test_split(
                df_dict[cell][0],
                test_size=validation_size_within,
                stratify=df_dict[cell][0].loc[:, stratify],
                random_state=random_state,
            )
        elif mtype == "across":
            df_train = df_dict["df"][df_dict["df"]["cell_type"] != cell]
            df_train, df_val = train_test_split(
                df_train,
                test_size=validation_size,
                stratify=df_train.loc[:, stratify],
                random_state=random_state,
            )

        X_train = df_train.drop(["label", "cell_type"], axis=1)
        y_train = df_train["label"]
        X_val = df_val.drop(["label", "cell_type"], axis=1)
        y_val = df_val["label"]
        del df_train
        del df_val

        best_params, fig = _optimize_parameters(
            [X_train, y_train, X_val, y_val],
            model_type,
            params,
            optim_time,
            n_trials,
            eval_metric,
            direction,
            random_state,
            cross_val=0,
        )

        params_opt_dict[cell] = best_params
        optim_fig_dict[cell] = fig

    return params_opt_dict, optim_fig_dict


def train_and_eval(
    df_dict: Dict[str, pd.DataFrame],
    mtype: str,
    model_type: str = "log_reg",
    params: dict = None,
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
    params = params or {}

    if not run:
        return (
            {"empty_model": ""},
            {"empty_model": plt.figure()},
            {"empty_model": pd.DataFrame()},
            {"empty_model": plt.figure()},
        )

    if run_name:
        mlflow.set_tag("mlflow.runName", run_name)

    start_time = time.time()
    shap_input = input(f"Do you want to generate SHAP plots for {model_type}? (y/n): ")
    if shap_input == "y":
        SHAP = True
    else:
        SHAP = False
    model_dict = _train_model(df_dict, mtype, model_type, params, SHAP=SHAP)
    metrics_dict_all, matrices_dict = _evaluate_model(
        model_dict, df_dict, mtype, model_type
    )
    feature_importances_dict, feature_importances_plot_dict = _get_feature_importances(
        model_dict, model_type
    )
    end_time = time.time()

    with open(f"data/08_reporting/Homo_sapiens/{mtype}/times.txt", "a") as file:
        file.write(
            f"Run_name: {run_name}\nTraining and evaluation of {model_type} for all cell types took: {str(round((end_time - start_time),2))} seconds.\n\n"
        )

    return (
        metrics_dict_all,
        matrices_dict,
        feature_importances_dict,
        feature_importances_plot_dict,
    )
