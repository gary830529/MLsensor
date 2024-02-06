# Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System

## Overview
This repository contains the R code accompanying the published paper "Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System." The paper and the code provide a comprehensive approach to predicting water quality using advanced machine learning techniques.

## Citation
If you use this code or our findings in your research, please cite:

Shyu, H.-Y., Castro, C. J., Bair, R. A., Lu, Q., & Yeh, D. H. (2023). Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System. ACS Environmental Au, 3(5), 308–318. https://doi.org/10.1021/acsenvironau.2c00072

## Abstract of Publication 
Developing advanced onsite wastewater treatment systems (OWTS) requires accurate and consistent water quality monitoring to evaluate treatment efficiency and ensure regulatory compliance. However, off-line parameters such as chemical oxygen demand (COD), total suspended solids (TSS), and *Escherichia coli* (*E. coli*) require sample collection and time-consuming laboratory analyses that do not provide real-time information of system performance or component failure. While real-time COD analyzers have emerged in recent years, they are not economically viable for onsite systems due to cost and chemical consumables. This study aimed to design and implement a real-time remote monitoring system for OWTS by developing several multi-input and single-output soft sensors. The soft sensor integrates data that can be obtained from well-established in-line sensors to accurately predict key water quality parameters, including COD, TSS, and *E. coli* concentrations. The temporal and spatial water quality data of an existing field-tested OWTS operated for almost two years (n = 56 data points) were used to evaluate the prediction performance of four machine learning algorithms. These algorithms, namely, partial least square regression (PLS), support vector regression (SVR), cubist regression (CUB), and quantile regression neural network (QRNN), were chosen as candidate algorithms for their prior application and effectiveness in wastewater treatment predictions. Water quality parameters that can be measured in-line, including turbidity, color, pH, NH<sub>4</sub><sup>+</sup>, NO<sub>3</sub><sup>–</sup>, and electrical conductivity, were selected as model inputs for predicting COD, TSS, and E. coli. The results revealed that the trained SVR model provided a statistically significant prediction for COD with a mean absolute percentage error (MAPE) of 14.5% and R<sup>2</sup> of 0.96. The CUB model provided the optimal predictive performance for TSS, with a MAPE of 24.8% and R<sup>2</sup> of 0.99. None of the models were able to achieve optimal prediction results for E. coli; however, the CUB model performed the best with a MAPE of 71.4% and R<sup>2</sup> of 0.22. Given the large fluctuation in the concentrations of COD, TSS, and *E. coli* within the OWTS wastewater dataset, the proposed soft sensor models adequately predicted COD and TSS, while *E. coli* prediction was comparatively less accurate and requires further improvement. These results indicate that although water quality datasets for the OWTS are relatively small, machine learning-based soft sensors can provide useful predictive estimates of off-line parameters and provide real-time monitoring capabilities that can be used to make adjustments to OWTS operations.

## Features
- Implementation of multiple machine learning algorithms
- Data preprocessing and analysis scripts
- Model evaluation and validation techniques
- Visualizations for model performance and data insightsFeatures
- Implementation of multiple machine learning algorithms
- Data preprocessing and analysis scripts
- Model evaluation and validation techniques
- Visualizations for model performance and data insights

## Data

The dataset used in this project is available upon request due to privacy or confidentiality considerations. If you are interested in accessing the data for research or educational purposes, please follow these steps:

1. Send an email to dhyeh@usf.edu with the subject line "Data Request for MLsensor."
2. In your email, please include your affiliation, the purpose of your data request, and how you plan to use the data.
3. We will review your request and respond with further instructions on how to access the data.

Please note that the data is subject to certain usage restrictions and may require agreeing to a data usage agreement.

## License
<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://github.com/gary830529/MLsensor">Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://gary830529.github.io/Homepage/">Hsiang-Yang Shyu</a> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">Attribution 4.0 International<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>



