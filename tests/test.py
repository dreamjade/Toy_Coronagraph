import toycoronagraph.main as toy

#load the example target image (example.fits) in Toy_Coronagraph/toycoronagraph/example_data/ folder
toy_target = toy.Target()

#plot the preimage, and the image will auto save in origin.png
toy_target.plot_origin()

#plot the final image through vortex coronagraph with charge = 2, and the image will auto save in charge2_final.png
charge = 2
toy_target.plot_final(charge)

#now add a planet
toy_target.add_planet(pos=[100,0], brightness=0.001)

#plot the final image with planets
#plot_final(self, charge, coronagraph_type='vortex', add_planet=True, img_pixel=512, psf_range=16, rot_number=360, plot_dpi=300)
charge = 2
toy_target.plot_final(charge)