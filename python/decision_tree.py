import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import tree
import matplotlib.pyplot as plt
import JSON

df = pd.read_json('data/data.json')  
# for testing of algorithm:
data = {
    'sequence_length': [100, 150, 120, 130, 140, 110],
    'methylation': [0.8, 0.3, 0.5, 0.7, 0.2, 0.4],
    'is_cancerous': [1, 0, 1, 1, 0, 0]  # 1 for cancerous, 0 for non-cancerous
}
df = pd.DataFrame(data)

X = df.drop('is_cancerous', axis=1)
y = df['is_cancerous']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Decision Tree classifier
# Using entropy as the criterion for information gain
clf = DecisionTreeClassifier(criterion='entropy', random_state=42)

# Train the model
clf.fit(X_train, y_train)

# Make predictions
y_pred = clf.predict(X_test)

# Evaluate the model
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("\nClassification Report:\n", classification_report(y_test, y_pred))

# Plot the tree
plt.figure(figsize=(12,8))
tree.plot_tree(clf, filled=True, feature_names=X.columns, class_names=['Not Cancerous', 'Cancerous'])
plt.savefig('plots/decision_tree.png')
