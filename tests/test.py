import toycoronagraph as toy

#load the example target image (example.fits) in Toy_Coronagraph/toycoronagraph/example_data/ folder
toy_target = toy.Target()

#plot the preimage, and the image will auto save in origin.png
toy_target.plot_origin()

#plot the final image through vortex coronagraph with charge = 2, and the image will auto save in charge2_final.png
charge = 2
toy_target.plot_final(charge)

