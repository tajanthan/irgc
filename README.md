# irgc
Iteratively Reweighted Graph Cut

-------------------------------------------------------------------------------------------
This code implements the IRGC and IRGC+expansion algorithms	described in the paper 

	"Iteratively Reweighted Graph Cut for Multi-label MRFs with Non-convex Priors", 
	Thalaiyasingam Ajanthan, Richard Hartley, Mathieu Salzmann and Hongdong Li,
	IEEE Conference on Computer Vision and Pattern Recognition,
	June 2015.

	*Code Assumptions:
		1. The MRF energy has the following form
			E(x) = \sum \theta_{i}(x_i) + \sum \theta_{ij} (x_i, x_j), 
			where \theta_{ij} (x_i, x_j) = \gamma_{ij} \theta(|x_i - x_j|).
		2. The code currently supports MRF with 4-connected grid structure only.
			nodes labelled from 0 --> width * height - 1, in a grid structure.
			E.g. width = 3, height = 2
			0 -- 1 -- 2
			|	 |	  |
			3 -- 4 -- 5

	This code is for research purposes only, if you want to use it for commercial purpose 
	please contact us.
	If you use this code, please consider citing the aforementioned paper 
	in any resulting publication.
	
	Contact: thalaiyasingam.ajanthan@nicta.com.au
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
Directory "alpha_exapnsion_west" contains the alpha-expansion code downloaded from
http://vision.csd.uwo.ca/code/gco-v3.0.zip

	We made some minor changes in the energy.h file to make it work with our 
	IRGC+expansion algorithm. Therefore we included it with our source code.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
Directory "maxflow" contains the max-flow code downloaded from
http://vision.csd.uwo.ca/code/maxflow-v3.01.zip

	We made some minor changes to reduce the memory requirement of the Ishikawa graph. 
	Therefore we included it with our source code.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
To assist the user, example.cpp, Makefile and sample data files are provided.

	Example usage 
	irgc.exe 10 10 5 <sample>/toy_unary_10_10_5.txt <sample>/toy_binary_4_10_10.txt TQ 2 1
		runs IRGC+expansion algorithm on 10x10 image with 5 labels, with 
		truncated quadratic pairwise potential - truncated at lambda = 2.
-------------------------------------------------------------------------------------------