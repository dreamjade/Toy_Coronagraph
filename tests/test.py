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


"""Test for SimpleClass. 
"""

import unittest


class  SimpleClass (object) : 
    """A simple class."""
    def __init__ (self, num1, num2) :         
        self.num1 = num1        
        self.num2 = num2

    def  sum (self) : 
        """Sum.
        Returns: 
            int, the sum of num1 and num2. 
        """ return self.num1 + self.num2
        ㄋ
    def  difference (self) : 
        """Difference.
        Returns: 
            int, the difference of num1 and num2. 
        """ if self.num1 < self.num2: return self.num2 - self.num1 return self.num1 - self.num2       

class  TestSimpleClass (unittest.TestCase) : 
    """Test of SimpleClass.""" 
    def setup (self) :         
        self.simple_class = SimpleClass( 1 , 2 )
        # def test_sum(self): 
        # """Test sum. 
        # """ 
        # self.assertEqual(self.simple_class.sum(), 3, 'test sum fail')    

    def  test_difference (self) : 
        """Test difference."""         
        self.assertEqual(self.simple_class.difference(), 1 , 'test difference fail' ) 
        # self.simple_class.num1, self.simple_class.num2 = \ 
        # self.simple_class.num2, self.simple_class.num1 
        # self.assertEqual(self.simple_class.difference(), 1, 
        # 'test difference fail')
        



if __name__ == '__main__' : 
    unittest.main()