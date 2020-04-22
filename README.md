# TUCANN: TUmor Characterization using Artificial Neural Networks

*Code repository for "Tumor Characterization using Unsupervised Learning of Mathematical Relations within Breast Cancer Data" by Cristian Axenie and Daria Kurz submitted at ENNS ICANN2020 (15th – 18th September 2020).*

**Codebase:**

datasets - the experimental datasets (csv files) and their source, each in separate directories
models   - codebase to run and reproduce the experiments

**Directory structure:**

model/.

* create_init_network.m       - init networks (SOM + HL)
* error_std.m                 - error std calculation function
* tumor_estimator_core.m      - main script to run the system
* model_rmse.m                - RMSE calculation function 
* model_sse.m                 - SSE calculation function
* parametrize_learning_law.m  - function to parametrize learning
* present_tuning_curves.m     - function to visualize SOM tuning curves
* randnum_gen.m               - weight initialization function
* tumor_growth_model_fit.m    - function implementing ODE models
* tumor_growth_models_eval.m  - main evaluation on runtime
* visualize_results.m         - visualize output and internals
* visualize_runtime.m         - visualize runtime



**Usage:**

* **model/tumor_estimator_core.m** - main function that runs the system and generates the runtime output file (mat file)
* **model/tumor_growth_models_eval.m** - evaluation and plotting function reading the runtime output file

