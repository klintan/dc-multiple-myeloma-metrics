## Dream Challenge Multiple Myeloma metrics implemented in Python

All metrics found in the metrics.r file, reimplemented in Python.  Many of them make use of Sklearn implemented metrics. 

Available metrics: 
- TimeROC
- Integrated AUC
- F1
- BAC
- MCC
- Weighted average across studies for each metric 


### Installation
```python setup.py install```

### Usage

```
import metrics

metrics.weightedAverage()
metrics.timeROC()
metrics.f1()
``` 

### Future improvements
- Adhere to the Sklearn-interface, to be able to evaluate your machine learning pipeline for your medical research in code.
- A bunch of the algorithms only support the default settings, extend to support more parameters

### Disclaimer
Still under development, feel free to fork and do a PR for any improvements and or extensions.
Beware that some of these could be incorrectly implemented and if you find anything please let me know or create a PR or github issue.

### License

MIT License, see LICENSE file 


