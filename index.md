# Iteratively Reweighted Graph Cut for Multi-label MRFs with Non-convex Priors

## Abstract
While widely acknowledged as highly effective in computer vision, multi-label MRFs with non-convex priors are difficult to optimize. To tackle this, we introduce an algorithm that iteratively approximates the original energy with an appropriately weighted surrogate energy that is easier to minimize. Our algorithm guarantees that the original energy decreases at each iteration. In particular, we consider the scenario where the global minimizer of the weighted surrogate energy can be obtained by a multi-label graph cut algorithm, and show that our algorithm then lets us handle of large variety of non-convex priors. We demonstrate the benefits of our method over state-of-the-art MRF energy minimization techniques on stereo and inpainting problems.

## Publication
Iteratively Reweighted Graph Cut for Multi-label MRFs with Non-convex Priors,
Thalaiyasingam Ajanthan, Richard Hartley, Mathieu Salzmann, and Hongdong Li, 
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), June 2016.

[pdf][1] [arxiv][2] [poster][3] [code][4]

[1]: /docs/irgc.pdf "pdf"
[2]: https://arxiv.org/abs/1411.6340 "arxiv"
[3]: /docs/irgc_poster.pdf "poster"
[4]: https://github.com/tajanthan/irgc "code"

