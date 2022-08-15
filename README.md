# Physionet-Challenge-Matlab-Solution

This repository contains Matlab solution for the Physionet Challenge - 2022.


Team Name: Melbourne Kangas.

Team Members : Zaria Imran, Ethan Grooby, Vinayaka Vivekananda Malgi, Chiranjibi Situala, Faezaeh Marzbanrad, Sunil Aryal.
 
MATLAB-specific instructions
You can use our MATLAB example classifier code (link) as a template. Consider cloning or downloading this repository, replacing our code with your code, and adding the updated files to your repository.


AUTHORS.txt, LICENSE.txt, README.md: Update as appropriate. Please include your authors. Unfortunately, our submission system will be unable to read your README file to change how we run your code.


train_model.m: Do not change this script. It calls your team_training_code.m script. We will not use the train_model.m script from your repository, so any change made to this code will not be included.


team_training_code.m: Update this script to create and save your model.


run_model.m: Do not change this script. It loads your model by calling load_model and runs your model by calling your team_testing_code function for each patient ID. We will not use the run_model.m script from your repository, so any change made to this code will not be included.


team_testing_code.m: Update this script to load and run your model weights and any parameters from files in your submission. It takes the data and your model as input and returns binary and probabilistic classifier outputs for each class as output.


Confirm that your code compiles and runs in MATLAB R2021b or R2022a (when available).


Push or upload your code to the root/base directory of the master branch of your repository.


We will download your code, compile it using the MATLAB compiler (mcc -m train_model.m -a . and mcc -m run_model.m -a .), and run it on our machines or Google Cloud.
