import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, explained_variance_score, max_error, mean_squared_log_error, median_absolute_error
import matplotlib.pyplot as plt

# Function to plot actual vs predicted values
def plot_actual_vs_predicted(y_true, y_pred, model_name):
    plt.figure(figsize=(8, 6))
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.plot(y_true, y_true, color='red', linestyle='--')  # Straight line representing perfect predictions
    plt.title(f"{model_name} - Actual vs Predicted")
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.grid(True)
    plt.show()
# Function to read data from Excel file
def read_data(file_path):
    return pd.read_excel(file_path)

# Function to perform train-test split
def perform_train_test_split(data, test_size=0.2, random_state=42):
    X = data.drop(columns=['denergy'])
    y = data['denergy']
    return train_test_split(X, y, test_size=test_size, random_state=random_state)

# Function to standardize data
def standardize_data(X_train, X_test):
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    return X_train_scaled, X_test_scaled

# Function to implement regression algorithms and calculate metrics
def implement_regression_algorithms(X_train, y_train, X_test, y_test):
    algorithms = {
        'Linear Regression': LinearRegression(),
        'Decision Tree': DecisionTreeRegressor(),
        'Random Forest': RandomForestRegressor()
    }
    results = {}
    for name, model in algorithms.items():
        model.fit(X_train, y_train)
        y_train_pred = model.predict(X_train)
        y_test_pred = model.predict(X_test)
        train_metrics = calculate_metrics(y_train, y_train_pred)
        test_metrics = calculate_metrics(y_test, y_test_pred)
        results[name] = {'Train Metrics': train_metrics, 'Test Metrics': test_metrics,'Actual Values': y_test,  # Store actual values
            'Predicted Values': y_test_pred  # Store predicted values
        }
        
    return results

# Function to calculate evaluation metrics
def calculate_metrics(y_true, y_pred):
    mse = mean_squared_error(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    evs = explained_variance_score(y_true, y_pred)
    me = max_error(y_true, y_pred)
    mae_median = median_absolute_error(y_true, y_pred)
    
    return {
        'MSE': mse,
        'MAE': mae,
        'R^2': r2,
        'Explained Variance Score': evs,
        'Max Error': me,
        'Median Absolute Error': mae_median
    }
# Modify print_results function to include plotting
def print_results(results):
    for name, metrics in results.items():
        print(f"{name}:")
        print("Train Metrics:")
        for metric, value in metrics['Train Metrics'].items():
            print(f"  {metric}: {value:.4f}")
        print("Test Metrics:")
        for metric, value in metrics['Test Metrics'].items():
            print(f"  {metric}: {value:.4f}")
        print()
        
        # Plot actual vs predicted values for Decision Tree and Random Forest
        if name in ['Decision Tree', 'Random Forest']:
             plot_actual_vs_predicted(metrics['Actual Values'], 
                                     metrics['Predicted Values'], name)


# Main function
def main():
    # Read data
    data = read_data("ex.xlsx")
    
    # Train-test split
    X_train, X_test, y_train, y_test = perform_train_test_split(data)
    
    # Standardize data
    X_train_scaled, X_test_scaled = standardize_data(X_train, X_test)
    
    # Implement regression algorithms
    results = implement_regression_algorithms(X_train_scaled, y_train, X_test_scaled, y_test)
    
    # Print results
    print_results(results)

if __name__ == "__main__":
    main()
