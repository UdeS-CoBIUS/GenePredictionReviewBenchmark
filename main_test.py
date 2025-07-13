from g3po_main.src.launch_prediction import launch_augustus, launch_helixer
from g3po_main.src.convert2gar import augustus_to_gar, helixer_to_gar
from g3po_main.src.compar_lvl_nuc import main_script
from g3po_main.src.compar_lvl_exon import main_script_exon
import os


class TestBenchmarkClass:
    
    def test_prediction_augustus(self):
        ret_val = launch_augustus(os.path.join(os.getcwd(), "g3po_main"), 150, "/usr/bin", True)
        assert ret_val == True
    
    def test_prediction_helixer(self):
        ret_val= launch_helixer(os.path.join(os.getcwd(), "g3po_main"), 150, ".", True)
        assert ret_val == True

    def test_conversion_predictions_augustus_to_gar(self):
        ret_val= augustus_to_gar(os.path.join(os.getcwd(), "g3po_main"), "augustus", 150, True)
        assert ret_val == True
        
    def test_conversion_predictions_helixer_to_gar(self):
        ret_val= helixer_to_gar(os.path.join(os.getcwd(), "g3po_main"), "helixer", 150, True)
        assert ret_val == True
        
    def test_comparison_nucleotide_level_for_augustus(self):
        ret_val= main_script(os.path.join(os.getcwd(), "g3po_main"), ["augustus"], [150], True)
        assert ret_val == True
        
    def test_comparison_nucleotide_level_for_helixer(self):
        ret_val= main_script(os.path.join(os.getcwd(), "g3po_main"), ["helixer"], [150], True)
        assert ret_val == True
        
    def test_comparison_exon_level_for_augustus(self):
        ret_val= main_script_exon(os.path.join(os.getcwd(), "g3po_main"), ["augustus"], True)
        assert ret_val == True
        
    def test_comparison_exon_level_for_helixer(self):
        ret_val= main_script_exon(os.path.join(os.getcwd(), "g3po_main"), ["helixer"], True)
        assert ret_val == True
        
    
    
        
    
        
    