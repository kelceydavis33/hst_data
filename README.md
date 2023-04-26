# hst_data

This repostiory tracks a project conducted under Dr. Kate Rubin at San Diego State University which I completed as an undergraduate student. We tackled the issue of abundant cosmic rays in photometry data from HST. This repository includes a .py file containing code that solved this issue and several .png image files showcasing code performance. Each .png file has a form that looks like this:

![image](https://user-images.githubusercontent.com/47015033/234673090-7fa02544-adbe-4078-a02d-67d02662c07e.png)


The leftmost panel is the final image after our code has run. The middle image is the original HST data, contaminated by cosmic rays. The rightmost panel shows the identified cosmic rays which were removed from the image.

Red circles indicate known stars, which we do not want to subtract from the images. This project ended exploring the question of how we could quantify what was a star and what was a cosmic ray. While we have some ideas of how to acomplish this, the project remains open as I have moved on to a new research group since. I beleive Dr. Rubin plans to let a future student explore this concept further.

Cosmic ray removal is machine learning based, primarily using the drizzlepac machine learning library.
