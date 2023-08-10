import pytest

from src.illumina_exon_corrector import IlluminaExonCorrector, VoidExonCorrector

class TestExonCorrection:
    
    @pytest.mark.parametrize("short_introns, long_exons",
                             [({(5154787, 5159231),(5159335, 5163585),(5163676, 5165837),(5163717, 5165837),(5165952, 5168251),(5168357, 5171292),(5171347, 5187612),(5187711, 5194499),(5194693, 5203306),(5203486, 5206034),(5206161, 5213972),(5214075, 5220170),
                             (5219146, 5324823),(5220285, 5232327)},
                             [(5154647, 5154786), (5159232, 5159334), (5163586, 5163675),(5165838, 5165951), (5168252, 5168356), (5171293, 5171346), (5187613, 5187710),(5194500, 5194692), (5203307, 5203485), (5206035, 5206160), (5213973, 5214074),
                             (5220171, 5220284), (5232328, 5232505)]),
                             ({(6298274, 6300182),(6300298, 6304184),(6304312, 6304452),(6304624, 6308485),(6308689, 6310175),(6310603, 6314328),(6314500, 6315018)},
                             [(6300221, 6300297), (6304185, 6304311), (6304453, 6304623), (6308486, 6308688),(6310176, 6310602), (6314329, 6314499), (6315019, 6315154)]),
                             (dict(),
                             [(10207833, 10209221), (10210692, 10210809), (10211784, 10211922), (10212801, 10212969), (10214943, 10214973)])
                             ])
    def test_no_correction(self, short_introns, long_exons):
        # not using short exons as those are never calculated
        # from the bam file I get introns directly
        corrector = IlluminaExonCorrector.from_data(short_introns)
        corrected_exons = corrector.correct_exons(long_exons)
        assert corrected_exons == long_exons
        
        
    @pytest.mark.parametrize("short_introns, long_exons, expected",
                             [({(11689307, 11851227),(11776574, 11820421),(11788720, 11851119),(11797197, 11850752),(11818316, 11868704),(11818325, 11868704),(11868791, 11900038),(11877909, 12147885),(11900116, 12003767),
                             (11900116, 12021186),(11926162, 12147885),(12003816, 12021186),(12021281, 12042269),(12042327, 12045146)},
                             [(11818191, 11818324), (11868705, 11868790), (11900039, 11900115), (12003768, 12003815), (12021187, 12021284), (12042270, 12042326), (12045147, 12045483)],
                             [(11818191, 11818324), (11868705, 11868790), (11900039, 11900115), (12003768, 12003815), (12021187, 12021280), (12042270, 12042326), (12045147, 12045483)]),
                             ({(15808064, 16021540),(15811004, 16002551),(15876182, 15877864),(15877961, 15883264),(15883387, 15889157),(15884847, 16136202),(15889308, 15889874),(15889988, 15901730),(15901791, 15903624),
                             (15903714, 15908405),(15908513, 15913144),(15910508, 16136202)},
                             [(15875937, 15876185), (15877865, 15877960), (15883265, 15883386), (15889158, 15889307), (15889879, 15889987), (15901731, 15901790), (15903625, 15903713), (15908406, 15908512), (15913145, 15913698)],
                             [(15875937, 15876181), (15877865, 15877960), (15883265, 15883386), (15889158, 15889307), (15889879, 15889987), (15901731, 15901790), (15903625, 15903713), (15908406, 15908512), (15913145, 15913698)]),
                             ({(16001998, 16237561),(16171584, 16171950),(16172161, 16172713),(16172824, 16172915),(16173054, 16173391)},
                             [(16171519, 16171583), (16171951, 16172160), (16172710, 16172823), (16172916, 16172987)],
                             [(16171519, 16171583), (16171951, 16172160), (16172714, 16172823), (16172916, 16172987)])
                             ])
    def test_only_offset(self, short_introns, long_exons, expected):
        corrector = IlluminaExonCorrector.from_data(short_introns)
        corrected_exons = corrector.correct_exons(long_exons)
        assert corrected_exons == expected
    
    
    @pytest.mark.parametrize("short_introns, long_exons, expected",
                             [({(4797249, 4898804),(4878206, 4878677),(4878710, 4898806),(4881433, 4900490),(4898873, 4900490),(4900539, 4902533),(4902592, 4911329),(4902605, 4907223),(4907298, 4909609),(4909712, 4911178),
                             (4909712, 5079593),(4911356, 4915185),(4911508, 5028541),(4911509, 5028533),(4911510, 5028533),(4911511, 5028533),(4911521, 5028533),(4915399, 5028533),(4915400, 5028533),(4915401, 5028531),(4915401, 5028533)},
                             [(4878046, 4878209), (4898782, 4898872), (4900491, 4900538), (4902534, 4902604), (4907224, 4907297), (4909610, 4909711), (4911179, 4911355), (4915186, 4916833)],
                             [(4878046, 4878205), (4878678, 4878709), (4898807, 4898872), (4900491, 4900538), (4902534, 4902604), (4907224, 4907297), (4909610, 4909711), (4911179, 4911355), (4915186, 4916833)]),
                             ({(11484862, 11582141),(11512061, 11678362),(11512062, 11678364),(11582290, 11588801),(11588934, 11615377),(11615507, 11633186),(11633216, 11658597),(11658771, 11666374),(11666542, 11818188)},
                             [(11484453, 11484861), (11582142, 11582289), (11588802, 11588933), (11615378, 11615528), (11658598, 11658770), (11666375, 11666541), (11667189, 11670635)],
                             [(11484453, 11484861), (11582142, 11582289), (11588802, 11588933), (11615378, 11615506), (11633187, 11633215), (11658598, 11658770), (11666375, 11666541), (11667189, 11670635)])
                             ])
    def test_only_skipped(self, short_introns, long_exons, expected):
        corrector = IlluminaExonCorrector.from_data(short_introns)
        corrected_exons = corrector.correct_exons(long_exons)
        assert corrected_exons == expected

class TestIntrons:
    
    @pytest.mark.parametrize("test_file, chromosome, start, end, expected_introns, expected_counts",
                             [("short_reads_toy/chr1.bam", "chr1", 11484453, 11670635, 
                             {(11484862, 11582141),(11512061, 11678362),(11512062, 11678364),(11582290, 11588801),(11588934, 11615377),(11615507, 11633186),(11633216, 11658597),(11658771, 11666374),(11666542, 11818188)},
                             {(11484861, 11582141): 489,(11512060, 11678362): 2,(11512061, 11678364): 6,(11582289, 11588801): 578,(11588933, 11615377): 646,(11615506, 11633186): 596,(11633215, 11658597): 576,
                             (11658770, 11666374): 598,(11666541, 11818188): 571}),
                             ("short_reads_toy/chr1.bam", "chr1", 4561154, 4565934, 
                             {(4313746, 4566513),(4404340, 4562619),(4431824, 4562619),(4513716, 4633993),(4533690, 4595836),(4562892, 4563322),(4562892, 4563994),(4562894, 4738843),(4563690, 4563994),(4563714, 4563994),
                             (4564087, 4565358),(4564087, 4566513),(4566033, 4788201)},
                             {(4313745, 4566513): 1,(4404339, 4562619): 1,(4431823, 4562619): 1,(4513715, 4633993): 1,(4533689, 4595836): 1,(4562891, 4563322): 841,(4562891, 4563994): 781,(4562893, 4738843): 1,
                             (4564086, 4565358): 1232,(4563689, 4563994): 214,(4563713, 4563994): 342,(4564086, 4566513): 161,(4566032, 4788201): 1})
                             ])
    def test_introns(self, test_file, chromosome, start, end, expected_introns, expected_counts):
        corrector = IlluminaExonCorrector(chromosome, start, end, test_file)
        assert corrector.short_introns == expected_introns
        assert corrector.counts == expected_counts
        
