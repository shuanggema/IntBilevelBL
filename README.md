Algorithms  
-------
compound_gbridge: Bi-level selection methods based on a bridge-type penalty and Lasso.

Maintainer
-------
Jin Liu   <jin.liu@duke-nus.edu.sg>


Publication
-------
Liu, J., Huang, J., Zhang, Y., Lan, Q., Rothman, N., Zheng, T., & Ma, S. (2014). Integrative analysis of prognosis data on multiple cancer subtypes. Biometrics, 70(3), 480-488.


Description
-------
For a specific gene, two levels of selection need to be conducted. The first is to determine whether it is associated with any subtype at all. This is achieved using a bridge-type penalty. For a gene associated with at least one subtype, the second level of selection is to determine which subtype(s) it is associated with. This is achieved using a Lasso type penalty. The composition of the two penalties can achieve the desired two-level selection.

Usage
-------
source("CallC_com_bridge.r")
com_bridge(x,y,group.inner,group.outer,lambda,gamma)