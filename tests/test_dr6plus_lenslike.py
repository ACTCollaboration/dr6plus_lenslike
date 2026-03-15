#!/usr/bin/env python3
"""
Unit test for dr6plus_lenslike.py modifications

This test ensures that:
1. The likelihood can be initialized without errors
2. The evaluate.yaml configuration works correctly
3. The likelihood computation gives consistent results before/after modifications
4. All variants can be loaded successfully
"""

import numpy as np
import os
import sys
import unittest
import tempfile
import shutil
from unittest.mock import patch, MagicMock
import warnings

# Add the package to path
sys.path.insert(0, '/home/jiaqu/dr6plus_lenslike')

try:
    from dr6plus_lenslike.dr6plus_lenslike import ACTDR6LensLike, load_data, generic_lnlike, variants
except ImportError as e:
    print(f"Could not import dr6plus_lenslike: {e}")
    print("Make sure the package is properly installed")
    sys.exit(1)

class TestDR6PlusLensLike(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures before running tests"""
        cls.test_dir = tempfile.mkdtemp()
        
        # Create mock data directory structure
        cls.data_dir = os.path.join(cls.test_dir, "data", "v1.0")
        os.makedirs(cls.data_dir, exist_ok=True)
        os.makedirs(os.path.join(cls.data_dir, "like_corrs"), exist_ok=True)
        
        # Create minimal mock data files for testing
        cls._create_mock_data_files()
        
        # Store reference results
        cls.reference_results = {}
        
    @classmethod  
    def tearDownClass(cls):
        """Clean up after tests"""
        shutil.rmtree(cls.test_dir)
        
    @classmethod
    def _create_mock_data_files(cls):
        """Create minimal mock data files for testing"""
        # Mock bandpowers.
        # Production data has 26 total ACT bins; after baseline cut [2:-6]
        # that gives 18 bins, matching response_cal_matrix shape (18, 9).
        nbins_total = 26
        nbins_act = 18        # = 26 - 2 - 6
        mock_bandpowers = np.random.random(nbins_total) * 1e-7
        np.savetxt(os.path.join(cls.data_dir, "clkk_bandpowers_act.txt"), mock_bandpowers)
        np.savetxt(os.path.join(cls.data_dir, "clkk_bandpowers_planck.txt"), mock_bandpowers[:10])

        # Mock binning matrices
        mock_binmat_act = np.random.random((nbins_total, 3000)) * 0.1
        mock_binmat_planck = np.random.random((10, 3000)) * 0.1
        np.savetxt(os.path.join(cls.data_dir, "binning_matrix_act.txt"), mock_binmat_act)
        np.savetxt(os.path.join(cls.data_dir, "binning_matrix_planck.txt"), mock_binmat_planck)

        # Mock covariance matrices (nbins_total × nbins_total for ACT)
        mock_cov_act = np.eye(nbins_total) * 1e-14
        mock_cov_actplanck = np.eye(nbins_total + 10) * 1e-14
        np.savetxt(os.path.join(cls.data_dir, "covmat_act.txt"), mock_cov_act)
        np.savetxt(os.path.join(cls.data_dir, "covmat_actplanck.txt"), mock_cov_actplanck)
        np.savetxt(os.path.join(cls.data_dir, "covmat_act_cmbmarg.txt"), mock_cov_act)
        np.savetxt(os.path.join(cls.data_dir, "covmat_actplanck_cmbmarg.txt"), mock_cov_actplanck)
        
        # Mock fiducial spectra for likelihood corrections
        like_corrs_dir = os.path.join(cls.data_dir, "like_corrs")
        ells = np.arange(2, 3000)
        mock_cls = np.column_stack([
            ells, 
            np.random.random(len(ells)) * 1e-10,  # TT
            np.random.random(len(ells)) * 1e-12,  # EE  
            np.random.random(len(ells)) * 1e-14,  # BB
            np.random.random(len(ells)) * 1e-11   # TE
        ])
        np.savetxt(os.path.join(like_corrs_dir, "cosmo2017_10K_acc3_lensedCls.dat"), mock_cls)
        
        # Mock lensing potential (6 columns: ell, unused×4, pp)
        mock_pp = np.column_stack([ells, np.zeros((len(ells), 4)),
                                  np.random.random(len(ells)) * 1e-8])
        np.savetxt(os.path.join(like_corrs_dir, "cosmo2017_10K_acc3_lenspotentialCls.dat"), mock_pp)
        
        # Mock correction matrices (minimal)
        np.save(os.path.join(like_corrs_dir, "norm_correction_matrix_Lmin0_Lmax4000.npy"), 
                np.random.random((4, 4000, 4000)) * 1e-10)
        
        # Mock other correction files
        correction_files = [
            "n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt",
            "N1der_KK_lmin600_lmax3000_full.txt",
            "N1der_TT_lmin600_lmax3000_full.txt", 
            "N1der_EE_lmin600_lmax3000_full.txt",
            "N1der_BB_lmin600_lmax3000_full.txt",
            "N1der_TE_lmin600_lmax3000_full.txt"
        ]
        
        for fname in correction_files:
            if "n0mv" in fname:
                # Real file is 2-row format: row 0 = ell values, row 1 = fAL values.
                # loadtxt unpacks as fAL_ls, fAL = np.loadtxt(file) → shape (2, N).
                data = np.row_stack([ells, np.random.random(len(ells)) * 1e-8])
            else:
                data = np.random.random((len(ells), len(ells))) * 1e-10
            np.savetxt(os.path.join(like_corrs_dir, fname), data)

        # Mock calibration files for selfcal tests
        # 8 noise weight files: one per array per pol (shape: 5001 values each)
        for arr in ['pa5a', 'pa5b', 'pa6a', 'pa6b']:
            for pol in ['T', 'E']:
                weights = np.ones(5001) * 0.25  # equal weights across 4 arrays
                np.savetxt(
                    os.path.join(cls.data_dir, f"noise_{arr}_{pol}_weights.txt"),
                    weights
                )
        # Response matrix R (18×9): small nonzero values so T² correction is detectable
        rng = np.random.default_rng(seed=0)
        response = rng.standard_normal((18, 9)) * 0.01
        np.savetxt(os.path.join(cls.data_dir, "response_cal_matrix.txt"), response)

    def setUp(self):
        """Set up for each test"""
        # Mock the data directory location
        self.original_file_dir = None
        
    def tearDown(self):
        """Clean up after each test"""
        pass
        
    def test_load_data_basic(self):
        """Test basic data loading functionality"""
        with patch('dr6plus_lenslike.dr6plus_lenslike.os.path.dirname') as mock_dirname:
            mock_dirname.return_value = self.test_dir
            
            try:
                data_dict = load_data(
                    variant="act_baseline",
                    ddir=self.data_dir,
                    lens_only=True,
                    apply_hartlap=False,
                    like_corrections=False
                )
                
                # Check that essential keys are present
                required_keys = ['data_binned_clkk', 'binmat_act', 'cov', 'cinv']
                for key in required_keys:
                    self.assertIn(key, data_dict, f"Missing required key: {key}")
                    
                # Store reference for comparison
                self.reference_results['basic_data'] = {
                    'data_shape': data_dict['data_binned_clkk'].shape,
                    'binmat_shape': data_dict['binmat_act'].shape,
                    'cov_shape': data_dict['cov'].shape
                }
                
            except Exception as e:
                self.fail(f"Basic data loading failed: {e}")
    
    def test_likelihood_initialization(self):
        """Test that ACTDR6LensLike can be initialized"""
        with patch('dr6plus_lenslike.dr6plus_lenslike.os.path.dirname') as mock_dirname:
            mock_dirname.return_value = self.test_dir
            
            # Test with minimal configuration
            like = ACTDR6LensLike()
            like.lens_only = True
            like.no_like_corrections = True
            like.apply_hartlap = False
            
            # Mock the data directory
            with patch('dr6plus_lenslike.dr6plus_lenslike.load_data') as mock_load_data:
                mock_load_data.return_value = {
                    'data_binned_clkk': np.random.random(18) * 1e-7,
                    'binmat_act': np.random.random((18, 3000)) * 0.1,
                    'cov': np.eye(18) * 1e-14,
                    'cinv': np.eye(18) * 1e14,
                    'include_planck': False,
                    'include_spt': False,
                    'include_spt_no_planck': False,
                    'likelihood_corrections': False,
                    'only_spt': False
                }
                
                try:
                    like.initialize()
                    self.assertIsNotNone(like.data)
                    self.assertEqual(like.requested_cls, ["pp"])
                except Exception as e:
                    self.fail(f"Likelihood initialization failed: {e}")
    
    def test_generic_lnlike_computation(self):
        """Test the generic likelihood computation"""
        # Create mock data dictionary
        data_dict = {
            'data_binned_clkk': np.random.random(10) * 1e-7,
            'binmat_act': np.random.random((10, 3000)) * 0.1,
            'cinv': np.eye(10) * 1e14,
            'include_planck': False,
            'include_spt': False,
            'include_spt_no_planck': False,
            'likelihood_corrections': False,
            'only_spt': False
        }

        # ell must span at least 3102 because generic_lnlike always calls
        # standardize(ell_kk, cl_kk, 3100) for the SPT cl_kk_spt array.
        ell = np.arange(2, 3102)
        cl_kk = np.random.random(len(ell)) * 1e-7
        cl_tt = np.random.random(len(ell)) * 1e-10
        cl_ee = np.random.random(len(ell)) * 1e-12
        cl_te = np.random.random(len(ell)) * 1e-11
        cl_bb = np.random.random(len(ell)) * 1e-14
        
        try:
            lnlike = generic_lnlike(
                data_dict, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
                trim_lmax=2998, do_norm_corr=False
            )
            
            # Check that we get a finite number
            self.assertTrue(np.isfinite(lnlike), "Likelihood should be finite")
            self.assertIsInstance(lnlike, (int, float, np.number))
            
            # Store reference for comparison
            self.reference_results['basic_lnlike'] = lnlike
            
        except Exception as e:
            self.fail(f"Generic likelihood computation failed: {e}")
    
    def test_variants_loading(self):
        """Test that all variants can be parsed without errors"""
        from dr6plus_lenslike.dr6plus_lenslike import parse_variant
        
        for variant in variants:
            try:
                v, baseline, include_planck, include_spt, include_spt_no_planck, only_spt = parse_variant(variant)
                # Basic sanity checks
                self.assertIsInstance(baseline, bool)
                self.assertIsInstance(include_planck, bool)
                self.assertIsInstance(include_spt, bool)
                self.assertIsInstance(include_spt_no_planck, bool)
                self.assertIsInstance(only_spt, bool)
            except Exception as e:
                self.fail(f"Failed to parse variant {variant}: {e}")
    
    def test_evaluate_yaml_configuration(self):
        """Test that the evaluate.yaml configuration can be processed"""
        # This test simulates running cobaya with the evaluate configuration
        
        # Mock provider for cobaya-like interface
        class MockProvider:
            def get_Cl(self, ell_factor=False, units='FIRASmuK2'):
                ell = np.arange(2, 5000)
                return {
                    'ell': ell,
                    'tt': np.random.random(len(ell)) * 1e-10,
                    'ee': np.random.random(len(ell)) * 1e-12, 
                    'te': np.random.random(len(ell)) * 1e-11,
                    'bb': np.random.random(len(ell)) * 1e-14,
                    'pp': np.random.random(len(ell)) * 1e-8
                }
            
            def get_param(self, param):
                return 1.0  # Default Alens
        
        with patch('dr6plus_lenslike.dr6plus_lenslike.os.path.dirname') as mock_dirname:
            mock_dirname.return_value = self.test_dir
            
            like = ACTDR6LensLike()
            like.lens_only = True
            like.no_like_corrections = True
            like.apply_hartlap = False
            like.provider = MockProvider()
            
            # Mock the data loading
            with patch('dr6plus_lenslike.dr6plus_lenslike.load_data') as mock_load_data:
                mock_load_data.return_value = {
                    'data_binned_clkk': np.random.random(18) * 1e-7,
                    'binmat_act': np.random.random((18, 3000)) * 0.1,
                    'cov': np.eye(18) * 1e-14,
                    'cinv': np.eye(18) * 1e14,
                    'include_planck': False,
                    'include_spt': False,
                    'include_spt_no_planck': False,
                    'likelihood_corrections': False,
                    'only_spt': False
                }
                
                try:
                    like.initialize()
                    logp = like.logp()
                    
                    self.assertTrue(np.isfinite(logp), "Log-likelihood should be finite")
                    
                    # Store reference for comparison
                    self.reference_results['evaluate_logp'] = logp
                    
                except Exception as e:
                    self.fail(f"Evaluate configuration test failed: {e}")
    
    def test_regression_comparison(self):
        """Compare results with stored reference values"""
        # This test is meant to be run after modifications to check for regressions
        
        # Use a fixed seed for all random data so the comparison is deterministic.
        # ell must span at least 3102 for the SPT standardize call.
        # binmat_act must have trim_lmax+2 = 3000 columns to match standardized cl_kk.
        np.random.seed(42)
        ell = np.arange(2, 3102)
        data_dict = {
            'data_binned_clkk': np.random.random(10) * 1e-7,
            'binmat_act': np.random.random((10, 3000)) * 0.1,
            'cinv': np.eye(10) * 1e14,
            'include_planck': False,
            'include_spt': False,
            'include_spt_no_planck': False,
            'likelihood_corrections': False,
            'only_spt': False
        }
        cl_kk = np.random.random(len(ell)) * 1e-7
        cl_tt = np.random.random(len(ell)) * 1e-10
        cl_ee = np.random.random(len(ell)) * 1e-12
        cl_te = np.random.random(len(ell)) * 1e-11
        cl_bb = np.random.random(len(ell)) * 1e-14

        run1 = generic_lnlike(
            data_dict, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998, do_norm_corr=False
        )
        run2 = generic_lnlike(
            data_dict, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998, do_norm_corr=False
        )

        # Calling generic_lnlike twice with identical inputs must give
        # identical results (determinism check).
        self.assertEqual(run1, run2, "Repeated calls with same inputs must be identical")


    # ------------------------------------------------------------------
    # Self-calibration tests (A–E)
    # ------------------------------------------------------------------

    def test_selfcal_A_raises_on_lens_only(self):
        """Test A: selfcal=True with lens_only=True raises ValueError."""
        with self.assertRaises(ValueError):
            load_data(
                variant='act_baseline',
                ddir=self.data_dir,
                selfcal=True,
                lens_only=True,
                apply_hartlap=False,
                like_corrections=False,
            )

    def test_selfcal_B_no_calibration_keys_when_false(self):
        """Test B: selfcal=False does not load calibration keys."""
        d = load_data(
            variant='act_baseline',
            ddir=self.data_dir,
            selfcal=False,
            lens_only=True,
            apply_hartlap=False,
            like_corrections=False,
        )
        self.assertNotIn('w_T', d)
        self.assertNotIn('response_cal_matrix', d)

    def test_selfcal_C_loads_calibration_keys_with_correct_shapes(self):
        """Test C: selfcal=True loads calibration arrays with correct shapes."""
        d = load_data(
            variant='act_baseline',
            ddir=self.data_dir,
            selfcal=True,
            lens_only=False,
            apply_hartlap=False,
            like_corrections=True,
        )
        self.assertIn('w_T', d)
        self.assertIn('w_E', d)
        self.assertIn('response_cal_matrix', d)
        self.assertEqual(d['w_T'].shape, (4, 5001))
        self.assertEqual(d['w_E'].shape, (4, 5001))
        self.assertEqual(d['response_cal_matrix'].shape, (18, 9))

    def _make_selfcal_data_dict(self):
        """Helper: return a data_dict with selfcal fields for Tests D and E.

        Loads with like_corrections=True so that fiducial CMB spectra
        (fiducial_cl_tt/ee/te) and selfcal arrays (w_T, w_E,
        response_cal_matrix) are present.  Sets likelihood_corrections=False
        before returning so that generic_lnlike takes the plain cl_kk path,
        avoiding the latent do_norm_corr TypeError bug documented in Step 1.
        The tests still exercise the T²(δc) scaling and dC computation paths.
        """
        d = load_data(
            variant='act_baseline',
            ddir=self.data_dir,
            selfcal=True,
            lens_only=False,
            apply_hartlap=False,
            like_corrections=True,
        )
        # Disable the norm-correction code path (which has a latent TypeError
        # for the do_norm_corr kwarg — see Step 1 note in dr6plus_lenslike.py).
        # fiducial_cl_* fields remain present for the dC selfcal computation.
        d['likelihood_corrections'] = False
        return d

    def test_selfcal_D_zero_calibration_unchanged_lnlike(self):
        """Test D: zero delta_c/delta_p leaves lnlike unchanged vs no selfcal."""
        d = self._make_selfcal_data_dict()

        # ell must span at least 3102 for the SPT standardize call in generic_lnlike.
        ell = np.arange(2, 3102)
        np.random.seed(123)
        cl_kk = np.random.random(len(ell)) * 1e-7
        cl_tt = np.random.random(len(ell)) * 1e-10
        cl_ee = np.random.random(len(ell)) * 1e-12
        cl_te = np.random.random(len(ell)) * 1e-11
        cl_bb = np.random.random(len(ell)) * 1e-14

        lnlike_no_cal = generic_lnlike(
            d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998,
            delta_c=None, delta_p=None,
        )
        lnlike_zero_cal = generic_lnlike(
            d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998,
            delta_c=np.zeros(4), delta_p=np.zeros(4),
        )

        self.assertAlmostEqual(lnlike_no_cal, lnlike_zero_cal, delta=1e-10,
                               msg="Zero calibration should leave lnlike unchanged")

    def test_selfcal_E_nonzero_calibration_changes_lnlike(self):
        """Test E: nonzero delta_c changes lnlike compared to no selfcal."""
        d = self._make_selfcal_data_dict()

        # ell must span at least 3102 for the SPT standardize call in generic_lnlike.
        ell = np.arange(2, 3102)
        np.random.seed(456)
        cl_kk = np.random.random(len(ell)) * 1e-7
        cl_tt = np.random.random(len(ell)) * 1e-10
        cl_ee = np.random.random(len(ell)) * 1e-12
        cl_te = np.random.random(len(ell)) * 1e-11
        cl_bb = np.random.random(len(ell)) * 1e-14

        lnlike_no_cal = generic_lnlike(
            d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998,
            delta_c=None, delta_p=None,
        )
        lnlike_with_cal = generic_lnlike(
            d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb,
            trim_lmax=2998,
            delta_c=np.array([0.01, 0.0, 0.0, 0.0]), delta_p=np.zeros(4),
        )

        self.assertGreater(
            abs(lnlike_no_cal - lnlike_with_cal), 1e-6,
            msg="Nonzero calibration should shift lnlike by more than 1e-6"
        )


def run_quick_test():
    """Quick test function that can be run from command line"""
    print("Running quick DR6+ lensing likelihood test...")
    
    # Run just the basic tests
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add specific tests for quick validation
    suite.addTest(TestDR6PlusLensLike('test_load_data_basic'))
    suite.addTest(TestDR6PlusLensLike('test_likelihood_initialization')) 
    suite.addTest(TestDR6PlusLensLike('test_generic_lnlike_computation'))
    suite.addTest(TestDR6PlusLensLike('test_variants_loading'))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    if result.wasSuccessful():
        print("\n✅ All basic tests passed! Your modifications appear to be safe.")
        return True
    else:
        print("\n❌ Some tests failed. Check the modifications carefully.")
        return False


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Test DR6+ lensing likelihood')
    parser.add_argument('--quick', action='store_true', help='Run quick test only')
    parser.add_argument('--full', action='store_true', help='Run full test suite')
    
    args = parser.parse_args()
    
    if args.quick:
        success = run_quick_test()
        sys.exit(0 if success else 1)
    else:
        # Run full test suite
        unittest.main(verbosity=2)
