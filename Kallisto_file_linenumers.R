#install.packages("renv")

renv::init(force=TRUE, bioconductor = TRUE)

#renv::install("data.table")

# Load necessary library
library(data.table)  # fread is faster than read.table, especially for large files

# Define the directory where the files are stored (change this path to your directory)
dir_path <- "/sbgenomics/project-files/"

# Get a list of all files in the directory that match the pattern
file_list <- list.files(dir_path, pattern = "kallisto.abundance.tsv.gz", full.names = TRUE)

length(file_list)

# Initialize an empty data frame to store the results
result_df <- data.frame(filename = character(), rows = numeric(), columns = numeric(), stringsAsFactors = FALSE)

# Loop through each file
for (file in file_list) {
  # Read the file using fread (which handles large files efficiently)
  data <- fread(file)
  
  # Get the number of rows and columns
  num_rows <- nrow(data)
  num_cols <- ncol(data)
  
  # Append the results to the result dataframe
  result_df <- rbind(result_df, data.frame(filename = basename(file), rows = num_rows, columns = num_cols))
}

# View the result
print(result_df)