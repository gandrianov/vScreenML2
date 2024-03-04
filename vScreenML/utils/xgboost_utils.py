import argparse
import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.metrics import matthews_corrcoef, roc_auc_score

def train_model():

    parser = argparse.ArgumentParser()
    parser.add_argument("-features", required=True)
    parser.add_argument("-output-prefix", required=True)
    parser.add_argument("-nsplits", default=5)

    args = parser.parse_args()

    features = args.features
    prefix = args.output_prefix
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

    with open(f"{prefix}_columns.csv", "w") as fwr:
        columns = data.select_dtypes(['number']).drop("Class", axis=1).columns
        columns = ",".join(list(columns))
        fwr.write(columns)

    model.save_model(f"{prefix}_model.json")
    
def predict_vscreenml_score():

    parser = argparse.ArgumentParser()
    parser.add_argument("-features", required=True)
    parser.add_argument("-output", required=True)
    
    args = parser.parse_args()

    columns = ["TotalExposedSasa", "TotalBSA", "InterfaceHydrophobicSasa", 
               "InterfacePolarSasa", "InteractionScore", "FaAtrInteraction", 
               "FaRepInteraction", "FaSolInteraction", "FaElecInteraction", 
               "HBondBbScInteraction", "HBondScInteraction", "GenBonded", 
               "HBInterface", "InterfaceUnsat", "NH1", "NH2", "NH3", "OH", 
               "carbonyl O", "carboxylate O", "SideFlexAlpha", "SideFlexBeta", 
               "SideFlexOther", "BackFlexAlpha", "BackFlexBeta", 
               "BackFlexOther", "PiPi", "TStacking", "CationPi", "SaltBridge", 
               "TotalElec", "TotalHBond", "TotalHphobics", "6.6", "6.7", "6.8", 
               "6.16", "7.6", "7.7", "7.8", "7.16", "8.6", "8.7", "8.8", "8.16",
               "9.6", "9.7", "9.8", "9.16", "15.6", "15.7", "15.8", "15.16",
               "16.6", "16.7", "16.8", "16.16", "17.6", "17.7", "17.8", "17.16",
               "35.6", "35.7", "35.8", "35.16", "53.6", "53.7", "53.8", "53.16",
               "exactmw", "amw", "lipinskiHBA", "lipinskiHBD",
               "NumRotatableBonds", "NumHBD", "NumHBA", "NumHeavyAtoms",
               "NumAtoms", "NumHeteroatoms", "NumAmideBonds", "FractionCSP3",
               "NumRings", "NumAromaticRings", "NumAliphaticRings",
               "NumSaturatedRings", "NumHeterocycles", 
               "NumAromaticHeterocycles", "NumSaturatedHeterocycles",
               "NumAliphaticHeterocycles", "NumSpiroAtoms", 
               "NumBridgeheadAtoms", "labuteASA", "tpsa", "CrippenClogP", 
               "CrippenMR", "chi0v", "chi1v", "chi2v", "chi3v", "chi4v", 
               "chi0n", "chi1n", "chi2n", "chi3n", "chi4n", "hallKierAlpha", 
               "kappa1", "kappa2", "kappa3", "Phi", "NumHAcceptors", 
               "NumHDonors", "MolLogP", "Proximal", 
               "Hydrogen bond", "Ionic", "Salt bridge", "Cation-pi", 
               "Hydrophobic", "Halogen bond", "Repulsive", 
               "Water-bridged hydrogen bond",  "Amide-aromatic stacking", 
               "Weak hydrogen bond", "Covalent bond", "Atom overlap", 
               "Van der Waals clash", "Van der Waals", "Chalcogen bond", 
               "Chalcogen-pi", "Halogen-pi", "Orthogonal multipolar", 
               "Parallel multipolar", "Antiparallel multipolar", 
               "Tilted multipolar", "Multipolar", "Cation-nucleophile", 
               "Anion-electrophile", "Unfavorable anion-nucleophile", 
               "Unfavorable cation-electrophile", 
               "Unfavorable nucleophile-nucleophile", 
               "Unfavorable electrophile-electrophile", "Pi-stacking", 
               "Face-to-face pi-stacking", "Face-to-edge pi-stacking", 
               "Face-to-slope pi-stacking", "Edge-to-edge pi-stacking", 
               "Edge-to-face pi-stacking", "Edge-to-slope pi-stacking", 
               "Displaced face-to-face pi-stacking", 
               "Displaced face-to-edge pi-stacking", 
               "Displaced face-to-slope pi-stacking", "C_RESIDUE", "INERTIA_3", 
               "SMALLEST_SIZE", "SURFACE_HULL", "VOLUME_HULL", 
               "hydrophobic_kyte", "hydrophobicity_pocket", "p_Ccoo", 
               "p_N_atom", "p_Ooh", "p_aliphatic_residue", "p_aromatic_residue", 
               "p_negative_residue"]

    important_cols = ['TotalBSA', 'InteractionScore', 'FaRepInteraction', 
                      'FaSolInteraction', 'FaElecInteraction', 'NH1', 'NH2', 
                      'NH3', 'OH', 'carbonyl O', 'carboxylate O', '7.8', '7.16', 
                      '8.8', '16.6', '17.16', '35.8', 'NumRotatableBonds', 
                      'NumHeavyAtoms', 'NumHeteroatoms', 
                      'NumAromaticHeterocycles', 'NumSaturatedHeterocycles', 
                      'tpsa', 'kappa1', 'Hydrogen bond', 'Repulsive', 
                      'Weak hydrogen bond', 'C_RESIDUE', 'SMALLEST_SIZE', 
                      'hydrophobic_kyte', 'p_aromatic_residue', 
                      'p_negative_residue']

    all_features_clf   = XGBClassifier(use_label_encoder=False)
    imprt_features_clf = XGBClassifier(use_label_encoder=False)

    root = __file__
    root = root[:root.rfind("/")]
    root = root[:root.rfind("/")]

    all_features_clf.load_model(f"{root}/models/all_feats_model.json")
    imprt_features_clf.load_model(f"{root}/models/imprt_feats_model.json")

    features = pd.read_csv(args.features)

    features["vScreenML (all feats)"]   = all_features_clf.predict_proba(features[columns])[:,1]
    features["vScreenML (imprt feats)"] = imprt_features_clf.predict_proba(features[important_cols])[:,1]
    
    features.to_csv(args.output, index=False)
