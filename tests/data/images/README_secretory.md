# secretory.obo
Figures showing a subset of GO terms related to 'secretory.'

  * [**'Secretory' GO terms without relationships loaded**](#secretory-go-terms-without-relationships-loaded)    
    % go_plot --obo=tests/data/secretory.obo -o tests/data/images/secretory_r0.png

  * [**'Secretory' GO terms with relationships loaded**](#secretory-go-terms-with-relationships-loaded)    
    % go_plot --obo=tests/data/secretory.obo -o tests/data/images/secretory_r1.png -r


## 'Secretory' GO terms without relationships loaded
% go_plot --obo=tests/data/secretory.obo -o tests/data/images/secretory_r0.png
![secretory without relationships](secretory_r0.png)

## 'Secretory' GO terms with relationships loaded
% go_plot --obo=tests/data/secretory.obo -o tests/data/images/secretory_r1.png -r

| Color   | relationship         |
|---------|----------------------|
|magenta  | part_of              |
|purple   | regulates            |
|red      | positively regulates |
|blue     | negatively regulates |

![secretory with relationships](secretory_r1.png)

Copyright (C) 2010-2018. Haibao Tang et al. All rights reserved.
