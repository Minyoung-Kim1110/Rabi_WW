# Base image with miniconda
FROM continuumio/miniconda3

# Set working directory
WORKDIR /opt/project

# Copy your environment file (optional)
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml

# Activate the environment by default
SHELL ["conda", "run", "-n", "WW_Sim", "/bin/bash", "-c"]

