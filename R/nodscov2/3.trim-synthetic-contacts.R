################################################################################
##
##                Compare synthetic and real contact data
##
################################################################################

## Working environment----------------------------------------------------------
# Libraries
library(tidyverse)
library(foreach)
library(doParallel)
library(Rcpp)
rm(list = ls())

# Helper functions
source("R/nodscov2/helper-functions.R")
sourceCpp("cpp/helper-functions.cpp")

## Load real and synthetic contact data-----------------------------------------
conditions = expand.grid(network = c("poincare", "herriot"), sim = 1:10, full = c(T,F))

registerDoParallel(4)
verif = foreach (r=1:nrow(conditions)) %dopar% {
  r = 2
    n = as.character(conditions[r,"network"])
    s = as.character(conditions[r,"sim"])
    full_network = as.logical(conditions[r, "full"])
    
    # List of individuals 
    admission = read.csv2(paste0("data/data-synthetic-graphs/input/admission_", n, ".csv")) %>%
      mutate(
        firstDate = as_datetime(firstDate, tz = "CET"),
        lastDate = as_datetime(lastDate, tz = "CET")
      )
    
    # Schedule of all individuals from the original dta 
    agenda = read.csv2(paste0("data/data-synthetic-graphs/input/agenda_", n, ".csv")) %>%
      mutate(
        firstDate = as_datetime(firstDate, tz = "CET"),
        lastDate = as_datetime(lastDate, tz = "CET")
      )
    
    # Synthetic data
    reconstructed_path = paste0("data/data-synthetic-graphs/biased/", n, "/")
    if (full_network) reconstructed_path = paste0("data/data-synthetic-graphs/full/", n, "/")
    
    all_f = list.files(reconstructed_path, pattern = "^matContact.*csv", full.names = T)
    f = all_f[grepl(paste0(s, "_oneday"), all_f)]
    
    interactions = read.csv2(f) %>%
      mutate(
        date_posix = as_datetime(date_posix, tz = "CET"),
        nSim = gsub("^.*Networks|_oneday.*$", "", f)
      ) %>%
      mutate(n = 1:n()) %>%
      nest(.by = c(nSim, n)) %>%
      mutate(data = map(data, trim_interactions_agenda, admission, agenda)) %>%
      unnest(data) %>%
      filter(length > 0, before_schedule == F) %>%
      select(from, to, date_posix, length)
    
    write.csv2(interactions, paste0(reconstructed_path, "truncated_data_", s, ".csv"), row.names = F)
}
