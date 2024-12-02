import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, accuracy_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import joblib

data = pd.read_csv('/Users/melaniesosa/Desktop/FALL 2024/POST-GENOMICS /GROUP PROJECT/Dataset Used for KNN/clinical_data_predictive_model.csv')

# Data Collection
print(f"Initial shape of data: {data.shape}")
print("Column-wise missing values:\n", data.isnull().sum())

# Pre-processing
# Remove columns with >20% missing values
threshold = 0.2
missing_ratio = data.isnull().mean()
data = data.loc[:, missing_ratio <= threshold]

# Remove rows with missing values
data = data.dropna()
print(f"Data shape after removing missing values: {data.shape}")
print("Remaining missing values:\n", data.isnull().sum())

# Descriptive Stats
print("Summary statistics of cleaned data:\n", data.describe())
target_column = 'Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code'
print("Target variable distribution:\n", data[target_column].value_counts())

# Features and Targets
X = data.drop(columns=[target_column])
y = data[target_column]

# Remove rare classes from the target variable
class_counts = y.value_counts()
rare_classes = class_counts[class_counts < 2].index
data = data[~data[target_column].isin(rare_classes)]

# Reassign X and y after filtering rare classes
X = data.drop(columns=[target_column])
y = data[target_column]

# Convert categorical variables to numeric codes
for col in X.select_dtypes(include=['object']).columns:
    X[col] = X[col].astype('category').cat.codes

# Remove low-variance features
low_variance_cols = X.columns[X.var(axis=0) < 0.01]
X = X.drop(columns=low_variance_cols)

# Remove highly correlated features (threshold = 0.9)
correlation_matrix = X.corr().abs()
upper_triangle = np.triu(np.ones(correlation_matrix.shape), k=1)
high_corr_indices = np.where(correlation_matrix > 0.9)
high_corr_features = [X.columns[col] for row, col in zip(*high_corr_indices) if row < col]
X = X.drop(columns=high_corr_features)

print(f"Final number of features: {X.shape[1]}")
print(f"Remaining features: {X.columns.tolist()}")

# Model Training and Evaluation
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42, stratify=y)
print(f"Training data shape: {X_train.shape}")
print(f"Test data shape: {X_test.shape}")

# ---- Graph Before Scaling ----
k_values = range(1, 21)
accuracies = []

for k in k_values:
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(X_train, y_train)
    y_pred = knn.predict(X_test)
    accuracies.append(accuracy_score(y_test, y_pred))

# Plot accuracy vs k (before scaling)
plt.plot(k_values, accuracies, marker='o')
plt.title('KNN Accuracy vs k (Before Scaling)')
plt.xlabel('Number of Neighbors (k)')
plt.ylabel('Accuracy')
plt.grid(True)
plt.xticks(ticks=k_values)
plt.show()

# ---- Scaling and Retraining ----
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Graph After Scaling
scaled_accuracies = []

for k in k_values:
    knn_scaled = KNeighborsClassifier(n_neighbors=k)
    knn_scaled.fit(X_train_scaled, y_train)
    y_pred_scaled = knn_scaled.predict(X_test_scaled)
    scaled_accuracies.append(accuracy_score(y_test, y_pred_scaled))

# Plot accuracy vs k (after scaling)
plt.plot(k_values, scaled_accuracies, marker='o')
plt.title('KNN Accuracy vs k (After Scaling)')
plt.xlabel('Number of Neighbors (k)')
plt.ylabel('Accuracy')
plt.grid(True)
plt.xticks(ticks=k_values)
plt.show()

# Evaluating k Values After Scaling
knn_scaled_k5 = KNeighborsClassifier(n_neighbors=5)
knn_scaled_k5.fit(X_train_scaled, y_train)
accuracy_scaled_k5 = accuracy_score(y_test, knn_scaled_k5.predict(X_test_scaled))

knn_scaled_k6 = KNeighborsClassifier(n_neighbors=6)
knn_scaled_k6.fit(X_train_scaled, y_train)
accuracy_scaled_k6 = accuracy_score(y_test, knn_scaled_k6.predict(X_test_scaled))

knn_scaled_k7 = KNeighborsClassifier(n_neighbors=7)
knn_scaled_k7.fit(X_train_scaled, y_train)
accuracy_scaled_k7 = accuracy_score(y_test, knn_scaled_k7.predict(X_test_scaled))

knn_scaled_k8 = KNeighborsClassifier(n_neighbors=8)
knn_scaled_k8.fit(X_train_scaled, y_train)
accuracy_scaled_k8 = accuracy_score(y_test, knn_scaled_k8.predict(X_test_scaled))

print(f"Accuracy with k=5: {accuracy_scaled_k5:.2f}")
print(f"Accuracy with k=6: {accuracy_scaled_k6:.2f}")
print(f"Accuracy with k=7: {accuracy_scaled_k7:.2f}")
print(f"Accuracy with k=8: {accuracy_scaled_k8:.2f}")

# Confusion Matrix for Final Model
cm = confusion_matrix(y_test, knn_scaled_k5.predict(X_test_scaled))
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=knn_scaled_k5.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.title("Confusion Matrix (k=5, Scaled Features)")
plt.show()

# Reporting and Deployment
print("Final KNN model details:")
print(f"Number of neighbors (k): {knn_scaled_k5.n_neighbors}")
print(f"Accuracy on test set: {accuracy_scaled_k5:.2f}")

joblib.dump(knn_scaled_k5, 'knn_model_scaled_k5.pkl')

with open('knn_features_scaled.txt', 'w') as f:
    f.write("\n".join(X.columns))
print("Model and feature list saved for deployment.")
