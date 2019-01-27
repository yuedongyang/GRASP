# TH-GRASS: Accurate Prediction of Genome-wide RNA Secondary Structure Based on Extreme Gradient Boosting

## Getting Start
These instructions will get you a copy of the project up and running on your local machine.

### Environment

TH-GRASS has been implemented in `Python3`.

### Requirements
Installing requirements:  
```
pip3 install -r requirements.txt
```

or avoiding problems in multiple Python environments:  

```
python3 -m pip install -r requirements.txt
```

## Training

### Command line:  
```
    usage: train.py [-h] [--input INPUT] [--training_mode TRAINING_MODE]
                    [--dataset DATASET] [--window_size WINDOW_SIZE]
                    [--n_jobs N_JOBS]

    optional arguments:
    -h, --help            show this help message and exit
    --input INPUT         The full path of the input file.
    --training_mode TRAINING_MODE
                            There are two training mode: single or global,
                            'single' is trained by single dataset and 'global' is
                            trained by all datasets.
    --dataset DATASET     This is for dataset, 'py' is for PARS-Yeast, 'ph' is
                            for PARS-Human, 'pdb' is for NMR/X-ray and 'global' is
                            for the three dataset mixed.
    --window_size WINDOW_SIZE
                            The window size when truncating RNA sequences.
    --n_jobs N_JOBS       Number of jobs to run in parallel.
```

### Start Training: 
The input data has been pre-processed from raw dataset. [Here](https://github.com/sysu-yanglab/TH-GRASS/tree/master/preprocessing "Data preprocess") is the data preprocessing methods.  

Example of training:  
* Using one of the datasets(PARS-Yeast, PARS-Human or SS_PDB) as training:
```
    python3 train.py --input='./data/PARS_yeast/py_encode_37.csv' --training_mode='single' --dataset='py' --window_size=37
    python3 train.py --input='./data/PARS_human/ph_encode_37.csv' --training_mode='single' --dataset='ph' --window_size=37
    python3 train.py --input='./data/NMR_X-Ray/pdb_encode_37.csv' --training_mode='single' --dataset='pdb' --window_size=37
```
That you can obtain three single models and save them in the directory: "./model/".  

* Using all of the datasets as training:
```
    python3 train.py --training_mode="global" --dataset="global" --window_size=37
```
That you can obtain the consensus model and it has been saved in the directory: "./model".

### Output:
* model:  
    The models file have been saved in the `'./model/'` directory.  
    If you follow the above commands to training, you will get some models named like "py_37.model", where 'py' is the name of training dataset, '37' is the window size you select.

* parameters:   
    Because training with `GridSearchCV` and `StratifiedKFold` function, the best parameter combinations would be selected from the parameter candidates and it would be saved in the directory: `'./parameters/'`.

## Evaluation

The performance were mainly evaluated by AUC.  

### Command line:  
```
    usage: evaluation.py [-h] [--mode MODE] --test_dataset TEST_DATASET
                        --train_model TRAIN_MODEL [--window_size WINDOW_SIZE]
                        [--output OUTPUT]

    optional arguments:
    -h, --help            show this help message and exit
    --mode MODE           The trained model is trained by two training modes:
                            single or global, 'single' is trained by single
                            dataset and 'global' is trained by all datasets.
    --test_dataset TEST_DATASET
                            This is for test-dataset. Choose from: 'py' is PARS-
                            Yeast 'ph' is PARS-Human, 'pdb' is NMR/X-ray.
    --train_model TRAIN_MODEL
                            This is for using which trained model. Choose from:
                            'py' is PARS-Yeast 'ph' is PARS-Human, 'pdb' is
                            NMR/X-ray and 'global' is Consensus-model.
    --window_size WINDOW_SIZE
                            The window size when truncating RNA sequences.
    --output OUTPUT         The directory for saving predictions on test-dataset.
```

### Usage:
* Single models  
``` 
    python3 evaluation.py --mode='single' --test_dataset='ph' --train_model='py' --window_size=37
```

* Consensus model
```
    python3 evaluation.py --mode='global' --train_model='global' --window_size=37
```

## Testing
Also, you can submit your own test RNA-sequences.

### Input file  
FASTA format, please ensure that your sequences only contain A, C, G, T, and U.

### Usage  

* Command line:  
```
    usage: test.py [-h] --input_fasta INPUT_FASTA --model_file MODEL_FILE
                --window_size WINDOW_SIZE --output_file OUTPUT_FILE

    optional arguments:
    -h, --help            show this help message and exit
    --input_fasta INPUT_FASTA
                            Input the test file of RNA sequences(FASTA format,
                            please ensure that your sequences only contain A, C,
                            G, T, and U)
    --model_file MODEL_FILE
                            The full path of the trained model.
    --window_size WINDOW_SIZE
                            The window size when truncating RNA sequences.
    --output_file OUTPUT_FILE
                            The full path of the output file.
```

* Run:  
```
    python3 test.py --input_fasta=InputFile --model_file='./model/PARS_yeast/py_37.model' --window_size=37 --output_file='./output/preds'
```

## Cite
If you find this work useful in your research, please consider citing the paper:  
**"Accurate Prediction of Genome-wide RNA Secondary Structure Based on eXtreme Gradient Boosting."**


## Contact
`yuedong.yang@gmail.com`, `raojh6@mail2.sysu.edu.cn`  or `keyaobin@mail2.sysu.edu.cn`