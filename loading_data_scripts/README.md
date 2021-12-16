# Loading iEEG data across different recording system formats


### BCI2000

WashU, Loma Linda, Albany, Mayo Clinic all use this format.

Data collected at different time points can sometimes require different methods. For newer data, the most robust method seems to be using the mex files in the `BCI2000mex` folder within this folder. To use this method:

1. Download the `BCI2000mex` folder to wherever you are loading the data
2. If on a mac, allow not verified developer code by running this in the terminal:        
 `sudo xattr -r -d com.apple.quarantine /path/to/folder/BCI2000mex/*`
3. Open a clean instance of matlab. Make sure you do NOT have fieldtrip in your path
4. In matlab, run the following:
```
addpath(genpath('/path/to/folder/BCI2000mex'))
path_data = '/path/to/data/0000.dat'
[signal,states,parameters] = load_bcidat(path_data);
```

