All the raw data, saved models, and predictions are git-ignored to minimize the size of the repository.

### Project setup

There are three key files.

- `competition_scratch.qmd`: I fit all my models here. You can see my tuning procedures and model specifications.
- `final_paper.qmd`: I wrote the paper here. I also wrote code for some graphs and tables, and determined the final model accuracies.
- `feature_engineering.R`: All my feature engineering functions live here. `do_feature_engineering()` is the main method. First run it on the training data with `fe_type = "training"` (line 33 in `competition_scratch.qmd`), then run with `fe_type = "testing"` on the testing data. This ensures there is no data leakage in the linear interpolation for the morphological embeddings.  
Each other function in this script creates some group of features and is described in `final_paper.qmd`.

To fit models, I employ the [`tidymodels` framework](https://www.tidymodels.org). Their documentation and vignettes are comprehensive and should cover all methods I used.

### Running the repository

First, clone the repository. Then, if using RStudio (recommended), open the `478 competition.Rproj` file. Download the `neuron-synapse-prediction` folder of information from the Kaggle competition site and put it in the root directory. (For instance, the relative path to the feature weights would be `neuron-synapse-prediction/feature_weights.csv`). Then, run through the code chunks in `competition_scratch.qmd` to fit and save models, which can take some time.

Throughout `competition_scratch.qmd` and to knit `final_paper.qmd`, you'll need to change the paths of the saved models to the new versions. Alternatively, email me at [prayaggordy\@gmail.com](mailto:prayaggordy@gmail.com) and I'll share a ZIP file of my fitted models.
