# Load the package
library(readxl)

# Specify the path to your .xls file
file_path <- "DataVulnerabilityAppendix.xls"

# Read the file
data <- read_excel(file_path)

# View the first few rows of the data
head(data)

