# hst_data

This repostiory tracks a project conducted under Dr. Kate Rubin at San Diego State University which I completed as an undergraduate student. 

I hoped to explore data taken with
the Hubble Space Telescope (HST) that had issues with severe cosmic ray presence. The images, of M82, are the
only images taken with HST of the galaxy’s winds using an MG II filter, meaning the observations
contain information of critical importance about the MG II transition. This particular
transition line will help us probe gas structures that may differ from what has been revealed with other
transmission lines. 

We tackled the issue of abundant cosmic rays in photometry data from HST through python based machine learning. This repository includes a .py file containing code that solved this issue and several .png image files showcasing code performance. Each .png file has a form that looks like this:

![image](https://user-images.githubusercontent.com/47015033/234673090-7fa02544-adbe-4078-a02d-67d02662c07e.png)


The leftmost panel is the final image after our code has run. The middle image is the original HST data, contaminated by cosmic rays. The rightmost panel shows the identified cosmic rays which were removed from the image.

Red circles indicate known stars, which we do not want to subtract from the images. They are indicated by the region file included here. This project ended exploring the question of how we could quantify what was a star and what was a cosmic ray. While we have some ideas of how to acomplish this, the project remains open as I have moved on to a new research group since. I beleive Dr. Rubin plans to let a future student explore this concept further.

Cosmic ray removal is machine learning based, primarily using the drizzlepac machine learning library.
