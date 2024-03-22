from io import BytesIO

import pytest
from file_process.exceptions import DelimiterError
from numpy import nan
import pandas as pd

from file_process.txt.txt_processor import TxtFileProcessor
from tests.test_file_process import get_remote_file_obj, TXT_INPUT_FILES_PATH


class TestCSVFileProcessor:
    original_data_path = f'{TXT_INPUT_FILES_PATH}/txt_example.txt'

    def _get_file_and_remote_file_obj(self, path: str):
        file = open(path, 'rb')
        file_obj = BytesIO(file.read())
        return file, file_obj

    def test_read_file(self):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        res = TxtFileProcessor(file_bytes_io)
        assert isinstance(res.txt_data, str)

    def test_get_preview(self):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        file_processor = TxtFileProcessor(file_bytes_io)
        _, _, _, text_file_preview = file_processor.get_preview()

        assert text_file_preview == """IQGIZBXAZYBADFENHMYLCWYTCOLEZQTEBIMAFMGGOBMHIBCBMPCPXCXHHZLWGNMAYCXLCVMOFVNBRFPIGYNWDMOCMVHNUITUHGZKZFQ
HONIYRGFSELKIWFCKUUQQTFORXFVCFHDVUASLZOTLYDCSDFVUOBIOFQCONFSQIFWOAFMYBLBAISKUTRNDZGNVLUFGREFOYUZFLRRIUMQYTNXKEGPDYSWSOEPDXBPPMNTHQQUV
NKXGOOKSCAMBZHICNBERBOSAVNWNHZKQRYIGSSYXVRIDKLTBUUCRBIETKQPPWH
IDZQWHUHWNRRDAMZGSKHDXYRAOLWFQOKQCFOZFUGAHGUFFTIWFHOGICAUMZHEFEXARPYIGHSQTOAKBCVAUKIOWHSXWLWHHRBEQMWWXEHKRVLIAFACKWNALZPLRIYLASOVDFHMWVIULQQRXRONRPOCRIATOFYVCHABWDBH
KVOHKZBVHQBCFGTQLEVXUHOZTCYQBGMCZFROTWLMMIMIBQEKYRLSOUUMUXPTWHUKDIAVRPKGAKLOVGCUGBROTSBKWTEQFYPMSUAVTILHUQFNDNEMYOYXNAKYMOOOCTYYNFCWRKZNURYWZGFMFRATXNTFVEQ
WCKWXNOGLETGINPIHOHVFLCWQFDFWIAPNHQMQQGEOBGBGRDOSPDVNTSZEXGQXULKVXFBDIMEVFPT
HYXKWZVDOTPTBKYTBAFZADFUAHGHECAKSCOIDCKCQXCDIWGZVTQV
TUKCQCELLPYYLQHHZLNRHSFMXOXNTSEHLTFLXVKYXPMRUISWODSTQCHRDBMEHIHHKNSVCBBYVBXHQXHQKAUWDSMHQDQCDKPBEPZALQPCDGDMULLRFZRAIKLAILHLEIDANDEWKAMXILUOQNYEMCXSQVHFDXXUDBZBFOKZPAXO
KDEDLUWLVNQEQITPASXHUSNNMDLAKRZMZPQVAWKYRSYVEGWBDX
FEZGFYHRUOZMTHSYWYCBRLEBQHPLHRHUCTVOGBUMIPTGHVRDTGTEHFCSTAYQXICYTFCWSHRNFFBAESFS
TTOFMSBNXMSQCOFLCIIWB
UZPGDORLDTYFUGHOWXPZOAXMVGKLGGMKPCGPYUDMGPVSUWQFTPSTYNEMVLCPIVEOPKFAPZLRWTPDCWEESGHRNEWMWMZTTZXZTRMDDMKHFGNKCOGETYBYNTFXYXSGSILPDLCWCEWZVFDIXCYGUHKYKZCHNASUICQPTPXZHIOYAZNCFGXTNDNRXRAAPWBZVB
HRBHHSVRNQPHUNHGIHKFEQZF
ECTAOQUHZNQLVDDPQWZLYBQYHSXOMIBIDCSSSPIKGUQCKAGGFOYXGZKPZTDZHEWRWEDYAOGWWUYGSTOLXOUTELRXPGVAWLYOFDPIAHDAWPAVBCWYHZXXWXDEOSNTAXCKAEZFBRPPWTUMURZOXPLRZZTHVRDVXUNRTUOVOAQ
FORMONWBVXWEZWKMB"""
        # assert obs_preview == [
        #     {"sepal_length": 5.1, "sepal_width": 3.5, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 4.9, "sepal_width": 3.0, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 4.7, "sepal_width": 3.2, "petal_length": 1.3, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 4.6, "sepal_width": 3.1, "petal_length": 1.5, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 5.0, "sepal_width": 3.6, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 5.4, "sepal_width": 3.9, "petal_length": 1.7, "petal_width": 0.4, "species": "setosa"},
        #     {"sepal_length": 4.6, "sepal_width": 3.4, "petal_length": 1.4, "petal_width": 0.3, "species": "setosa"},
        #     {"sepal_length": 5.0, "sepal_width": 3.4, "petal_length": 1.5, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 4.4, "sepal_width": 2.9, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
        #     {"sepal_length": 4.9, "sepal_width": 3.1, "petal_length": 1.5, "petal_width": 0.1, "species": "setosa"}
        # ]
        # assert var_preview == None
        # assert var_names == ["sepal_length", "sepal_width", "petal_length", "petal_width", "species"]
