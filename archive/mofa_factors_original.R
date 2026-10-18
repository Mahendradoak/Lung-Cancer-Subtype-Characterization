# =============================================================================
# R SCRIPT TO GENERATE THE MOFA FACTORS FILE FOR PYTHON
# =============================================================================

# Set your working directory to the project folder
setwd("C:/Users/mahen/Desktop/research internship UCD/final project")

# Load the MOFA2 library
library(MOFA2)

# Load your pre-trained MOFA model
model_filename <- "trained_mofa_model.rds"

if (!file.exists(model_filename)) {
  stop("CRITICAL ERROR: The file 'trained_mofa_model.rds' was not found.")
}
MOFAobject <- readRDS(model_filename)

# Extract the latent factor scores
factor_scores <- get_factors(MOFAobject, factors = "all")[[1]]

# Save the factors to the CSV file
write.csv(factor_scores, "mofa_factors.csv", row.names = TRUE)

# Print a confirmation message
cat("Successfully created 'mofa_factors.csv'. You can now run your full Python script.")
