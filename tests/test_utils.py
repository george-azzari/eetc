"""
Untilities for unittests
"""

def compare_bands(testcase, img, expected_bands, assert_equal_kwargs=None):
    if assert_equal_kwargs is None:
        assert_equal_kwargs = {}

    img_bands = img.bandNames().getInfo()
    expected_bands = sorted(list(set(expected_bands)))
    img_bands = sorted(list(img_bands))

    testcase.assetEqual(img_bands, expected_bands, **assert_equal_kwargs)
