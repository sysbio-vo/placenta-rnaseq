This directory holds code and scripts used in differential expression analysis.

1. Get data by running `./get_data.sh` script
	1.1. Internaly it's goind to download raw data from gdrive with gdown.
	1.2. Copy in to raw_counts folder
	1.3. Copy merged counts data from aligment folder, if data is not avalible locally downloads from project's google dirve
	1.4. Merge data together, output to counts

