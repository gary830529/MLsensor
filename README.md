# Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System

## Overview
This repository contains the R code accompanying the published paper "Development of a Soft Sensor Using Machine Learning Algorithms for Predicting the Water Quality of an Onsite Wastewater Treatment System." The paper and the code provide a comprehensive approach to predicting water quality using advanced machine learning techniques.

## Abstract of Publication 
Developing advanced onsite wastewater treatment systems (OWTS) requires accurate and consistent water quality monitoring to evaluate treatment efficiency and ensure regulatory compliance. However, off-line parameters such as chemical oxygen demand (COD), total suspended solids (TSS), and Escherichia coli (E. coli) require sample collection and time-consuming laboratory analyses that do not provide real-time information of system performance or component failure. While real-time COD analyzers have emerged in recent years, they are not economically viable for onsite systems due to cost and chemical consumables. This study aimed to design and implement a real-time remote monitoring system for OWTS by developing several multi-input and single-output soft sensors. The soft sensor integrates data that can be obtained from well-established in-line sensors to accurately predict key water quality parameters, including COD, TSS, and E. coli concentrations. The temporal and spatial water quality data of an existing field-tested OWTS operated for almost two years (n = 56 data points) were used to evaluate the prediction performance of four machine learning algorithms. These algorithms, namely, partial least square regression (PLS), support vector regression (SVR), cubist regression (CUB), and quantile regression neural network (QRNN), were chosen as candidate algorithms for their prior application and effectiveness in wastewater treatment predictions. Water quality parameters that can be measured in-line, including turbidity, color, pH, NH4+, NO3–, and electrical conductivity, were selected as model inputs for predicting COD, TSS, and E. coli. The results revealed that the trained SVR model provided a statistically significant prediction for COD with a mean absolute percentage error (MAPE) of 14.5% and R2 of 0.96. The CUB model provided the optimal predictive performance for TSS, with a MAPE of 24.8% and R<sup>2</sup> of 0.99. None of the models were able to achieve optimal prediction results for E. coli; however, the CUB model performed the best with a MAPE of 71.4% and R2 of 0.22. Given the large fluctuation in the concentrations of COD, TSS, and E. coli within the OWTS wastewater dataset, the proposed soft sensor models adequately predicted COD and TSS, while E. coli prediction was comparatively less accurate and requires further improvement. These results indicate that although water quality datasets for the OWTS are relatively small, machine learning-based soft sensors can provide useful predictive estimates of off-line parameters and provide real-time monitoring capabilities that can be used to make adjustments to OWTS operations.