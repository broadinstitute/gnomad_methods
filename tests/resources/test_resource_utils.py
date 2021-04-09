from unittest.mock import patch

from gnomad.resources import resource_utils


class TestTableResource:
    @patch("hail.read_table")
    def test_read_table(self, read_table):
        resource = resource_utils.TableResource("gs://my-bucket/table.ht")

        ds = resource.ht()
        read_table.assert_called_with("gs://my-bucket/table.ht")
        assert ds == read_table.return_value


class TestMatrixTableResource:
    @patch("hail.read_matrix_table")
    def test_read_matrix_table(self, read_matrix_table):
        resource = resource_utils.MatrixTableResource("gs://my-bucket/matrix_table.mt")

        ds = resource.mt()
        read_matrix_table.assert_called_with("gs://my-bucket/matrix_table.mt")
        assert ds == read_matrix_table.return_value


class TestPedigreeResource:
    @patch("hail.Pedigree.read")
    def test_read_pedigree(self, read_pedigree):
        resource = resource_utils.PedigreeResource("gs://my-bucket/pedigree.ped")

        ds = resource.pedigree()
        read_pedigree.assert_called()
        print(read_pedigree.call_args)
        assert read_pedigree.call_args[0][0] == "gs://my-bucket/pedigree.ped"
        assert ds == read_pedigree.return_value

    @patch("hail.import_fam")
    def test_read_fam(self, import_fam):
        resource = resource_utils.PedigreeResource("gs://my-bucket/pedigree.fam")

        ds = resource.ht()
        import_fam.assert_called()
        assert import_fam.call_args[0][0] == "gs://my-bucket/pedigree.fam"
        assert ds == import_fam.return_value


class TestBlockMatrixResource:
    @patch("hail.linalg.BlockMatrix.read")
    def test_read_block_matrix(self, read_block_matrix):
        resource = resource_utils.BlockMatrixResource("gs://my-bucket/block_matrix.bm")

        ds = resource.bm()
        read_block_matrix.assert_called_with("gs://my-bucket/block_matrix.bm")
        assert ds == read_block_matrix.return_value
