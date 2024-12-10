#CSV file
df <- read.csv("~/Genomics/Final_Project/clinical_data_predictive_model.csv")
#________________REMOVE NA's___________________
# given there are 388 observations in total in some cases is better to remove the entire column
# if 80% of the data is missing the column will be removed
df = df[,!sapply(df, function(x) mean(is.na(x)))>0.2]
df = na.omit(df)
#check NA's
colSums(is.na(df)) #left with total of 206 obs and 44 variables
#________________FACTORIZE_____________________#
df[sapply(df, is.character)] = lapply(df[sapply(df, is.character)], as.factor)
summary(df)
#______________TRAINING DATA_______________________#
n = nrow(df) 
n_training = round(n * 0.5)
n_test = n - n_training
index = sample(n, size = n_training)
training_data = df[index ,] #50% of dataset
#________________TEST DATA____________________________#
test_data = df[-index ,] #50% of dataset
#______________Correlation Matrix_____________________#
library(dplyr)
cor_d <- cor(df[sapply(df, is.numeric)]) # Select only numeric columns
library(corrplot)
corrplot(cor_d, method = "color", tl.cex = 0.5, # Adjust text size
         col = colorRampPalette(c("blue", "white", "red"))(200))
threshold <- 0.3
cor_matrix_filtered <- cor_d
# Filter: Set correlations below the threshold to 0 (or NA, which would hide them in `corrplot`)
cor_matrix_filtered[abs(cor_matrix_filtered) < threshold] <- NA
# Visualize with corrplot
corrplot(cor_matrix_filtered, method = "color", na.label = " ", tl.cex = 0.5,
         col = colorRampPalette(c("blue", "white", "red"))(200))
#_________________RANDOM FOREST__________________________#
library(randomForest)
summary(test_data)
model_rf = randomForest(df$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code ~ ., data = df, ntree = 100)
#print model
summary(model_rf)
print(model_rf)
#importance features
importance(model_rf)
plot(model_rf)
#produce variable importance plot
varImpPlot(model_rf) 
#Metric Evaluation
predictions <- predict(model_rf, df)
library(caret)
confusionMatrix(predictions, df$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)
#___________________SVM MODEL_________________________
sapply(df, function(x) if (is.factor(x)) nlevels(x) else NA)

library(e1071)

# Remove columns with only one unique level
training_data <- training_data[sapply(training_data, function(col) length(unique(col)) > 1)]
# SVM formula
training_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code <- as.factor(training_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)

# Re-run the SVM
svmfit <- svm(Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code ~ ., 
              data = training_data, 
              kernel = "radial", 
              cost = 1, 
              gamma = 0.1, 
              scale = TRUE)
print(svmfit)
summary(svmfit)
#Prediction and Evaluation
predictions <- predict(svmfit, test_data)
confusion <- confusionMatrix(predictions, test_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)
print(confusion)

library(caret)
confusionMatrix(predictions, test_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)
#Feature Scaling
training_data_scaled <- scale(training_data[sapply(training_data, is.numeric)])

#________________KNN________________________________

#___________________MLP MODEL_________________________#
library(neuralnet)

# Ensure the target variable is numeric
training_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code <- 
  as.numeric(as.factor(training_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)) - 1

# Scale the predictors
training_data_scaled <- as.data.frame(lapply(training_data[sapply(training_data, is.numeric)], scale))
training_data_scaled$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code <- 
  training_data$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code

# Define the formula dynamically
formula <- as.formula(paste("Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code ~", 
                            paste(names(training_data_scaled)[-length(names(training_data_scaled))], 
                                  collapse = " + ")))

# Train the neural network
model <- neuralnet(formula,
                   data = training_data_scaled,
                   hidden = c(5, 3),
                   linear.output = FALSE)

# Plot the model
plot(model)

# Ensure test data has the same columns as training data
test_data_scaled <- test_data[, colnames(training_data_scaled), drop = FALSE]

# Remove the target variable from test data
test_data_scaled <- test_data_scaled[, -which(names(test_data_scaled) == "Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code")]

# Predict on the test data
predictions <- compute(model, test_data_scaled)

# Extract predicted values
predicted_classes <- ifelse(predictions$net.result > 0.5, 1, 0)                          



# Convert to factors and align levels
actual_classes <- as.factor(test_data_scaled$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)
predicted_classes <- factor(predicted_classes, levels = levels(actual_classes))

# Remove unused levels
actual_classes <- droplevels(actual_classes)
predicted_classes <- droplevels(predicted_classes)

# Ensure there are overlapping levels
if (length(intersect(levels(actual_classes), levels(predicted_classes))) == 0) {
  stop("No overlapping levels between predicted and actual classes. Check your data!")
}

# Compute confusion matrix
library(caret)
confusionMatrix(predicted_classes, actual_classes)
plot(svmfit)
