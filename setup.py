import setuptools

setuptools.setup(
    name="carps",
    version="0.1",
    author="Randy Munoz",
    author_email="randy.munoz@yale.edu",
    description="Tool for Calculating the Angle of Ram Pressure Stripping",
    packages=setuptools.find_packages(include=['Data_Analysis','CARPS_site']),
    python_requires='>=3',
    install_requires=["numpy","matplotlib", "streamlit", "astropy", "streamlit_image_coordinates", "astroquery","pandas"]
    )