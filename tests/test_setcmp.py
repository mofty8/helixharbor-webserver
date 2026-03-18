import unittest
from unittest.mock import patch

import pandas as pd

import app as helix_app
import function_set


def _tmh_df(protein_id, full_sequence="MAVVVVVVVVVVVVVVVVVV"):
    return pd.DataFrame(
        [
            {
                "id": protein_id,
                "type": "transmembrane region",
                "begin": 1,
                "end": 5,
                "full_sequence": full_sequence,
                "sequence": full_sequence[:5],
                "SA": 10.0,
                "volume": 20.0,
                "bulkiness": 30.0,
                "hydrophobicity": 1.5,
            }
        ]
    )


def _aac_df(protein_id, sequence="MAVVV"):
    return pd.DataFrame(
        [
            {
                "id": protein_id,
                "type": "transmembrane region",
                "begin": 1,
                "end": len(sequence),
                "full_sequence": sequence,
                "sequence": sequence,
                "SA": 10.0,
                "volume": 20.0,
                "bulkiness": 30.0,
                "hydrophobicity": 1.5,
            }
        ]
    )


class RunForListCmpTests(unittest.TestCase):
    def test_backfills_missing_full_sequence_for_both_lists(self):
        background = pd.DataFrame(columns=["id", "begin", "end", "full_sequence", "type"])

        def fake_retrieve(uniprot_id):
            return (
                [
                    {
                        "id": uniprot_id,
                        "type": "transmembrane region",
                        "begin": 1,
                        "end": 4,
                        "full_sequence": None,
                    }
                ],
                [],
            )

        with patch.object(function_set, "retrieve_protein_features", side_effect=fake_retrieve), patch.object(
            function_set,
            "fetch_uniprot_sequence",
            side_effect=lambda uniprot_id: "MAVVV",
        ):
            list1_df, list2_df = function_set.run_for_list_cmp(["P11111"], ["Q22222"], background)

        self.assertEqual(list1_df.iloc[0]["full_sequence"], "MAVVV")
        self.assertEqual(list2_df.iloc[0]["full_sequence"], "MAVVV")
        self.assertEqual(list1_df.iloc[0]["sequence"], "MAVV")
        self.assertEqual(list2_df.iloc[0]["sequence"], "MAVV")

    def test_mixed_background_and_uniprot_rows_do_not_create_duplicate_type_columns(self):
        background = pd.DataFrame(
            [
                {
                    "id": "P11111",
                    "Type": "transmembrane region",
                    "begin": 1,
                    "end": 4,
                    "full_sequence": "MAVVV",
                    "sequence": "MAVV",
                    "hydrophobicity": 1.2,
                }
            ]
        )

        def fake_retrieve(uniprot_id):
            return (
                [
                    {
                        "id": uniprot_id,
                        "type": "transmembrane region",
                        "begin": 1,
                        "end": 4,
                        "full_sequence": "VVVVM",
                    }
                ],
                [],
            )

        with patch.object(function_set, "retrieve_protein_features", side_effect=fake_retrieve):
            list1_df, list2_df = function_set.run_for_list_cmp(["P11111"], ["Q22222"], background)

        self.assertEqual(list(list1_df.columns).count("type"), 1)
        self.assertEqual(list(list2_df.columns).count("type"), 1)
        extracted = function_set.extract_second_row_values(list1_df, "option4", "list1")
        self.assertEqual(len(extracted), 1)


class SetComparisonRouteTests(unittest.TestCase):
    def setUp(self):
        helix_app.app.config["TESTING"] = True
        self.client = helix_app.app.test_client()

    def test_returns_helpful_error_when_lists_have_no_comparable_tmh_values(self):
        empty_compare_df = pd.DataFrame(columns=["id", "type", "begin", "end", "full_sequence", "sequence"])

        with patch.object(helix_app, "load_tmh_dbs", return_value=pd.DataFrame()), patch.object(
            helix_app, "run_for_list_cmp", return_value=(empty_compare_df, empty_compare_df)
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "setcmp",
                    "list_1": "P11111",
                    "list_2": "Q22222",
                    "option2": "option1",
                    "option3": "option4",
                    "sequence": "",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"No comparable transmembrane regions were found", response.data)

    def test_returns_plot_for_valid_compare_request(self):
        with patch.object(helix_app, "load_tmh_dbs", return_value=pd.DataFrame()), patch.object(
            helix_app,
            "run_for_list_cmp",
            return_value=(_tmh_df("P11111"), _tmh_df("Q22222", full_sequence="VVVVMMMMM")),
        ), patch.object(helix_app, "save_density_raw_data"), patch.object(
            helix_app, "denisty_plot", return_value="ZmFrZS1wbmc="
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "setcmp",
                    "list_1": "P11111",
                    "list_2": "Q22222",
                    "option2": "option1",
                    "option3": "option4",
                    "sequence": "",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,ZmFrZS1wbmc=", response.data)

    def test_returns_plot_for_all_tmh_compare_request(self):
        with patch.object(helix_app, "load_tmh_dbs", return_value=pd.DataFrame()), patch.object(
            helix_app,
            "run_for_list_cmp",
            return_value=(_tmh_df("P11111"), _tmh_df("Q22222", full_sequence="VVVVMMMMM")),
        ), patch.object(helix_app, "save_density_raw_data"), patch.object(
            helix_app, "denisty_plot", return_value="YWxsLXRtaC1wbmc="
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "setcmp",
                    "list_1": "P11111",
                    "list_2": "Q22222",
                    "option2": "option2",
                    "option3": "option4",
                    "sequence": "",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,YWxsLXRtaC1wbmc=", response.data)

    def test_returns_aac_plot_for_compare_request(self):
        with patch.object(helix_app, "load_tmh_dbs", return_value=pd.DataFrame()), patch.object(
            helix_app,
            "run_for_list_cmp",
            return_value=(_aac_df("P11111"), _aac_df("Q22222", sequence="VVVVM")),
        ), patch.object(helix_app, "aac_density_plot", return_value="YWFjLXBsb3Q="):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "setcmp",
                    "list_1": "P11111",
                    "list_2": "Q22222",
                    "option2": "option1",
                    "option3": "option6",
                    "sequence": "",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,YWFjLXBsb3Q=", response.data)

    def test_requires_both_lists(self):
        response = self.client.post(
            "/HelixHarbor/",
            data={
                "input_type": "setcmp",
                "list_1": "P11111",
                "list_2": "",
                "option2": "option1",
                "option3": "option4",
                "sequence": "",
            },
        )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Please provide both List 1 and List 2", response.data)


class ListAgainstBackgroundRouteTests(unittest.TestCase):
    def setUp(self):
        helix_app.app.config["TESTING"] = True
        self.client = helix_app.app.test_client()

    def test_returns_plot_for_first_tmh_background_comparison(self):
        query_df = _tmh_df("P11111")
        background_df = _tmh_df("BG1111")

        with patch.object(helix_app, "load_tmh_dbs", return_value=background_df.copy()), patch.object(
            helix_app, "run_for_tmh_list", return_value=(query_df.copy(), [])
        ), patch.object(helix_app, "save_density_raw_data"), patch.object(
            helix_app, "denisty_plot", return_value="YmctcGxvdA=="
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "uniprot",
                    "sequence": "P11111",
                    "option4": "option1",
                    "option": "option1",
                    "option2": "option1",
                    "option3": "option4",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,YmctcGxvdA==", response.data)

    def test_returns_plot_for_all_tmh_background_comparison(self):
        query_df = _tmh_df("P11111")
        background_df = _tmh_df("BG1111")

        with patch.object(helix_app, "load_tmh_dbs", return_value=background_df.copy()), patch.object(
            helix_app, "run_for_tmh_list", return_value=(query_df.copy(), [])
        ), patch.object(helix_app, "save_density_raw_data"), patch.object(
            helix_app, "denisty_plot", return_value="YWxsLWJnLXBsb3Q="
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "uniprot",
                    "sequence": "P11111",
                    "option4": "option1",
                    "option": "option1",
                    "option2": "option2",
                    "option3": "option4",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,YWxsLWJnLXBsb3Q=", response.data)

    def test_returns_aac_plot_for_background_comparison(self):
        query_df = _aac_df("P11111")
        background_df = _aac_df("BG1111", sequence="VVVVM")

        with patch.object(helix_app, "load_tmh_dbs", return_value=background_df.copy()), patch.object(
            helix_app, "run_for_tmh_list", return_value=(query_df.copy(), [])
        ), patch.object(helix_app, "aac_density_plot", return_value="YmctYWFj"):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "uniprot",
                    "sequence": "P11111",
                    "option4": "option1",
                    "option": "option1",
                    "option2": "option1",
                    "option3": "option6",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"data:image/png;base64,YmctYWFj", response.data)

    def test_returns_non_transmembrane_notice_with_plot(self):
        query_df = _tmh_df("P11111")
        background_df = _tmh_df("BG1111")

        with patch.object(helix_app, "load_tmh_dbs", return_value=background_df.copy()), patch.object(
            helix_app, "run_for_tmh_list", return_value=(query_df.copy(), ["Q99999"])
        ), patch.object(helix_app, "save_density_raw_data"), patch.object(
            helix_app, "denisty_plot", return_value="bm90aWNlLXBsb3Q="
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "uniprot",
                    "sequence": "P11111",
                    "option4": "option1",
                    "option": "option1",
                    "option2": "option1",
                    "option3": "option4",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Q99999", response.data)
        self.assertIn(b"data:image/png;base64,bm90aWNlLXBsb3Q=", response.data)

    def test_returns_helpful_error_when_background_comparison_has_no_tmh_values(self):
        empty_compare_df = pd.DataFrame(columns=["id", "type", "begin", "end", "full_sequence", "sequence"])

        with patch.object(helix_app, "load_tmh_dbs", return_value=empty_compare_df.copy()), patch.object(
            helix_app, "run_for_tmh_list", return_value=(empty_compare_df.copy(), [])
        ):
            response = self.client.post(
                "/HelixHarbor/",
                data={
                    "input_type": "uniprot",
                    "sequence": "P11111",
                    "option4": "option1",
                    "option": "option1",
                    "option2": "option1",
                    "option3": "option4",
                },
            )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"No comparable transmembrane regions were found for your list or the selected background", response.data)


if __name__ == "__main__":
    unittest.main()
