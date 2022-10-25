import argparse
import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.metrics import matthews_corrcoef, roc_auc_score

def train_model():

    parser = argparse.ArgumentParser()
    parser.add_argument("-features", required=True)
    parser.add_argument("-model-name", required=True)
    parser.add_argument("-nsplits", default=5)

    args = parser.parse_args()

    features = args.features
    model_name = args.model_name
    n_splits = int(args.nsplits)

    data = pd.read_csv(features)
    data = data.sample(frac=1).reset_index(drop=True)
    data = data.select_dtypes(['number'])

    X = data.drop("Class", axis=1).values
    y = data["Class"]

    model = XGBClassifier(use_label_encoder=False)
    skf = StratifiedKFold(n_splits=n_splits)

    prob_predictions = []

    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        model.fit(X_train, y_train, eval_metric="logloss")
        
        df = pd.DataFrame(data.iloc[test_index,:])
        df["Proba"] = model.predict_proba(X_test)[:,1]
        
        prob_predictions.append(df)

    prob_predictions = pd.concat(prob_predictions)

    acc = accuracy_score(prob_predictions["Class"], round(prob_predictions["Proba"],0))
    prc = precision_score(prob_predictions["Class"], round(prob_predictions["Proba"],0))
    rec = recall_score(prob_predictions["Class"], round(prob_predictions["Proba"],0))
    mcc = matthews_corrcoef(prob_predictions["Class"], round(prob_predictions["Proba"],0))
    auc_score = roc_auc_score(prob_predictions["Class"], prob_predictions["Proba"])

    print(f"Accuracy: {round(acc, 2)}")
    print(f"Precision: {round(prc, 2)}")
    print(f"Recall: {round(rec, 2)}")
    print(f"MCC: {round(mcc, 2)}")
    print(f"ROC AUC: {round(auc_score, 2)}")

    model.fit(X, y, eval_metric="logloss")

    root = __file__
    root = root[:root.rfind("/")]
    root = root[:root.rfind("/")]

    with open(f"{root}/models/{model_name}_columns.csv", "w") as fwr:
        columns = data.select_dtypes(['number']).drop("Class", axis=1).columns
        columns = ",".join(list(columns))
        fwr.write(columns)

    model.save_model(f"{root}/models/{model_name}.json")


def predict_vscreenml_score():

    parser = argparse.ArgumentParser()
    parser.add_argument("-features", required=True)
    parser.add_argument("-output", required=True)
    parser.add_argument("-model-name", default="DUDE")

    args = parser.parse_args()

    features_filename = args.features
    model_name = args.model_name
    output = args.output

    features = pd.read_csv(features_filename)

    root = __file__
    root = root[:root.rfind("/")]
    root = root[:root.rfind("/")]

    model = XGBClassifier(use_label_encoder=False)
    model.load_model(f"{root}/models/{model_name}.json")

    columns = open(f"{root}/models/{model_name}_columns.csv", "r").read().split(",")

    for c in columns:
        if c not in features.columns:
            raise Exception(f"Feature {c} is not presented in the {features_filename} file")

    features["Predicted_Class"] = model.predict(features[c])
    features["VScreenML_Score"] = model.predict_proba(features[c])[:,1]

    features.to_csv(output, index=False)