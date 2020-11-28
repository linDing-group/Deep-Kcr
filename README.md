# Deep-Kcr
Identification of Kcr sites using deep learning
accurate detection of lysine crotonylation sites using deep learning method. If you don't have enough computer background to run the code we posted, you can send your query sequences to my email (haolv@std.uestc.edu.cn) at any time.


The ‘code’ file contains all of the code for feature extraction and model training and testing.
1. ‘featureForCKSAAP_PWAA_AAindex_CTD_EBGW.py’ is used to extract features of CKSAAP, PWAA, AAindex, CTD, and EBGW. 
Usage: before using it, users need to prepare the temp file in advance and set the file name in the corresponding position in the code.
And input ‘python featureForCKSAAP_PWAA_AAindex_CTD_EBGW.py’ in cmd box.
2. ‘word2vec.ipynb’ is used to extract embedding kmer features.
3. ‘CNN_model.ipynb’ is used to train the convolutional neural network-based model.
4. ‘load_model.ipynb’ is used to get testing result on independent dataset.


Due to file size limitation of Github, all temp files can be downloaded from http://lin-group.cn/server/Deep-Kcr/download.html.
