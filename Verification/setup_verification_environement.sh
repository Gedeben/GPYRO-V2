#!/bin/bash

# Step 1: Create ~/envs directory if it doesn't exist
mkdir -p ~/envs

# Step 2: Create the Python virtual environment only if it doesn't exist
ENV_PATH=~/envs/env_gpyro_VandV
if [ ! -d "$ENV_PATH" ]; then
    python3 -m venv "$ENV_PATH"
    ENV_CREATED=true
else
    ENV_CREATED=false
fi

# Step 3: Activate the environment
source "$ENV_PATH/bin/activate"

# Step 4: Installrequirements

pip3 install -r requirements.txt

# Final message
if [ "$ENV_CREATED" = true ]; then
    echo "Environment 'env_gpyro_VandV' created and activated"
else
    echo "env_gpyro_VandV activated"
fi

