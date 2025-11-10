# Base image with Python and R
FROM rocker/tidyverse:4.3.2

# Install system deps
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev \
    libxml2-dev libcurl4-openssl-dev libssl-dev \
    && apt-get clean

# Set working directory
WORKDIR /app

# Copy everything (scripts, data, configs)
COPY . /app/


# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Install R dependencies
RUN Rscript install_R_packages.R

# Make scripts executable
RUN find /app -name "*.py" -exec chmod +x {} \;

# Set default command 
CMD ["bash"]
