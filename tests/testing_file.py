from glidesign.geometry import Airfoil
    
if __name__ == '__main__':
    from parapy.gui import display
    wing = Airfoil(
        airfoil_name = 'NACA 24012',
        chord = 1.5,
        cst_poly_order = 5,
        num_points = 200
    )
    display(wing)