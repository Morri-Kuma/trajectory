"""
tests/test_lineage_graphsim_sctimebench.py

Unit tests for lineage_graphsim_sctimebench.py asserting scTimeBench parity.

Coverage:
  1.  modified_floyd_warshall -- max-min bottleneck reachability
  2.  floyd_warshall_closure  -- binary transitive closure with diagonal=1
  3.  Jaccard uses full-matrix flatten (no diagonal exclusion)
  4.  Confusion metrics use full-matrix flatten (no diagonal exclusion)
  5.  all_paths closure diagonal=1 in both pred and ref (TP on diagonal)
  6.  all_paths AUC uses original W_pred (diagonal zeroed), not W_reach
  7.  Threshold selection via sklearn PRC/ROC on full flatten
  8.  simple criterion end-to-end
  9.  all_paths criterion end-to-end
  10. compute_graph_sim_metrics schema
  11. lineage_metrics.json schema from eval_lineage.run_lineage_evaluation
  12. align_predicted_to_reference
  13. AUC diagonal zeroing (scTimeBench parity)
  14. all_paths thresholds original W_pred, not W_reach

Run with:
    cd C:\\Users\\37620\\trajectory
    python -m pytest benchmark/evaluation/tests/test_lineage_graphsim_sctimebench.py -v
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from benchmark.evaluation.lineage_graphsim_sctimebench import (
    modified_floyd_warshall,
    floyd_warshall_closure,
    align_predicted_to_reference,
    compute_graph_sim_metrics,
    validate_predicted_label_space,
    _select_threshold,
    _compute_simple_metrics,
    _compute_all_paths_metrics,
    _jaccard_flat,
    _confusion_metrics,
)


# ---------------------------------------------------------------------------
# 1. modified_floyd_warshall
# ---------------------------------------------------------------------------

class TestModifiedFloydWarshall:
    def test_no_edges(self):
        W = np.zeros((4, 4))
        assert np.allclose(modified_floyd_warshall(W), 0.0)

    def test_direct_edge_preserved(self):
        W = np.zeros((3, 3))
        W[0, 1] = 0.8; W[1, 2] = 0.5
        R = modified_floyd_warshall(W)
        assert R[0, 1] == pytest.approx(0.8)
        assert R[1, 2] == pytest.approx(0.5)

    def test_indirect_better_than_direct(self):
        # 0->2 direct=0.3, 0->1=0.9, 1->2=0.7 => bottleneck via 1 = min(0.9,0.7)=0.7 > 0.3
        W = np.zeros((3, 3))
        W[0, 2] = 0.3; W[0, 1] = 0.9; W[1, 2] = 0.7
        assert modified_floyd_warshall(W)[0, 2] == pytest.approx(0.7)

    def test_chain_bottleneck(self):
        W = np.zeros((4, 4))
        W[0, 1] = 0.6; W[1, 2] = 0.4; W[2, 3] = 0.9
        assert modified_floyd_warshall(W)[0, 3] == pytest.approx(0.4)

    def test_no_backward_path(self):
        W = np.zeros((3, 3))
        W[0, 1] = 0.7; W[1, 2] = 0.5
        assert modified_floyd_warshall(W)[2, 0] == pytest.approx(0.0)

    def test_output_shape(self):
        assert modified_floyd_warshall(np.eye(5) * 0.5).shape == (5, 5)


# ---------------------------------------------------------------------------
# 2. floyd_warshall_closure -- scTimeBench diagonal=1 semantics
# ---------------------------------------------------------------------------

class TestFloydWarshallClosure:
    def test_diagonal_always_one(self):
        B = np.zeros((4, 4), dtype=int)
        B[0, 1] = 1
        C = floyd_warshall_closure(B)
        assert np.all(np.diag(C) == 1), "Diagonal must be 1 (dist[i][i]=0 < inf in FW)"

    def test_isolated_node_diagonal_still_one(self):
        B = np.zeros((3, 3), dtype=int)
        assert np.all(np.diag(floyd_warshall_closure(B)) == 1)

    def test_transitive_chain(self):
        B = np.zeros((4, 4), dtype=int)
        B[0, 1] = B[1, 2] = B[2, 3] = 1
        C = floyd_warshall_closure(B)
        for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
            assert C[i, j] == 1, f"Expected reachability C[{i},{j}]=1"
        assert C[3, 0] == 0

    def test_disconnected_components(self):
        B = np.zeros((4, 4), dtype=int)
        B[0, 1] = 1; B[2, 3] = 1
        C = floyd_warshall_closure(B)
        assert C[0, 2] == 0 and C[0, 3] == 0
        assert C[0, 0] == 1 and C[2, 2] == 1

    def test_output_dtype_int(self):
        B = np.array([[0, 1], [0, 0]], dtype=int)
        assert floyd_warshall_closure(B).dtype in (np.int32, np.int64, int)


# ---------------------------------------------------------------------------
# 3. Jaccard uses full-matrix flatten (no diagonal exclusion)
# ---------------------------------------------------------------------------

class TestJaccardFullMatrix:
    def test_perfect_match_includes_diagonal(self):
        n = 4
        ref = np.zeros((n, n), dtype=int)
        ref[0, 1] = ref[1, 2] = 1
        np.fill_diagonal(ref, 1)
        pred = ref.copy()
        assert _jaccard_flat(pred.flatten(), ref.flatten()) == pytest.approx(1.0)

    def test_empty_union_returns_nan(self):
        assert math.isnan(_jaccard_flat(np.zeros(4), np.zeros(4)))

    def test_partial_overlap(self):
        pred = np.array([1, 1, 0, 0])
        ref  = np.array([1, 0, 1, 0])
        assert _jaccard_flat(pred, ref) == pytest.approx(1.0 / 3.0)


# ---------------------------------------------------------------------------
# 4. Confusion metrics use full-matrix flatten
# ---------------------------------------------------------------------------

class TestConfusionMetricsFullMatrix:
    def test_diagonal_ones_count_as_tp(self):
        n = 3
        ref  = np.ones((n, n), dtype=int)
        pred = np.ones((n, n), dtype=int)
        cm = _confusion_metrics(pred.flatten(), ref.flatten())
        assert cm["precision"] == pytest.approx(1.0)
        assert cm["recall"]    == pytest.approx(1.0)
        assert cm["f1"]        == pytest.approx(1.0)

    def test_perfect_prediction(self):
        ref = np.array([0, 1, 0, 0, 0, 1, 0, 0, 0])
        cm  = _confusion_metrics(ref, ref)
        assert cm["precision"] == pytest.approx(1.0)
        assert cm["recall"]    == pytest.approx(1.0)

    def test_all_wrong(self):
        ref  = np.array([1, 1, 0, 0])
        pred = np.array([0, 0, 1, 1])
        cm   = _confusion_metrics(pred, ref)
        assert cm["precision"] == pytest.approx(0.0)
        assert cm["recall"]    == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# 5. all_paths closure diagonal=1 in both pred and ref (TP on diagonal)
# ---------------------------------------------------------------------------

class TestAllPathsDiagonalParity:
    def _chain(self):
        n = 4
        B = np.zeros((n, n), dtype=int)
        B[0, 1] = B[1, 2] = B[2, 3] = 1
        W = np.full((n, n), 0.02)
        W[0, 1] = W[1, 2] = W[2, 3] = 0.8
        np.fill_diagonal(W, 0.0)
        return W, B

    def test_ref_closure_diagonal_is_one(self):
        _, B = self._chain()
        assert np.all(np.diag(floyd_warshall_closure(B)) == 1)

    def test_pred_closure_diagonal_is_one(self):
        W, _ = self._chain()
        W_reach = modified_floyd_warshall(W)
        pred_direct = (W_reach >= 0.3).astype(int)
        assert np.all(np.diag(floyd_warshall_closure(pred_direct)) == 1)

    def test_diagonal_tp_in_all_paths_metrics(self):
        """
        With 0.02 noise background, PRC selects threshold=0.02 so all entries
        are predicted, giving recall=1.0. Jaccard = B_reach.sum() / n^2 = 10/16.
        This is correct scTimeBench behavior: PRC uses full-flatten labels before
        FW closure and cannot anticipate that closure recovers diagonal TP.
        """
        W, B = self._chain()
        W_reach = modified_floyd_warshall(W)
        B_reach = floyd_warshall_closure(B)
        bundle = _compute_all_paths_metrics(
            W, W_reach, B_reach, 0.1, True, True
        )
        assert bundle["recall"] == pytest.approx(1.0)
        n = 4
        expected_jaccard = float(B_reach.sum()) / (n * n)
        assert bundle["jaccard_similarity"] == pytest.approx(expected_jaccard, abs=1e-9)


# ---------------------------------------------------------------------------
# 6. all_paths AUC uses original W_pred (with diagonal zeroed), not W_reach
# ---------------------------------------------------------------------------

class TestAllPathsAUCUsesOriginalPred:
    def test_auc_from_w_pred_not_w_reach(self):
        n = 4
        W_pred = np.zeros((n, n))
        W_pred[0, 1] = 0.9; W_pred[1, 2] = 0.7; W_pred[2, 3] = 0.5
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = B_ref[2, 3] = 1

        W_reach = modified_floyd_warshall(W_pred)
        B_reach = floyd_warshall_closure(B_ref)

        bundle = _compute_all_paths_metrics(W_pred, W_reach, B_reach, 0.1, True, True)

        from sklearn.metrics import roc_auc_score
        y_true  = B_reach.flatten().astype(int)
        # Expected AUC: W_pred with diagonal zeroed (already 0 here)
        W_pred_dz = W_pred.copy(); np.fill_diagonal(W_pred_dz, 0.0)
        if len(np.unique(y_true)) >= 2:
            expected = float(roc_auc_score(y_true, W_pred_dz.flatten()))
            assert bundle["auc_roc"] == pytest.approx(expected, abs=1e-9)

    def test_w_reach_differs_from_w_pred_for_indirect(self):
        n = 3
        W_pred = np.array([[0.0, 0.9, 0.0], [0.0, 0.0, 0.8], [0.0, 0.0, 0.0]])
        W_reach = modified_floyd_warshall(W_pred)
        assert W_reach[0, 2] == pytest.approx(0.8)
        assert W_pred[0, 2]  == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# 7. Threshold selection via sklearn PRC/ROC on full flatten
# ---------------------------------------------------------------------------

class TestThresholdSelectionSklearn:
    def test_prc_selects_correct_threshold(self):
        n = 3
        ref   = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
        pred  = np.full((n, n), 0.05)
        pred[0, 1] = pred[1, 2] = 0.8
        t, tt = _select_threshold(pred.flatten(), ref.flatten(), prc_threshold=True)
        assert tt == "prc"
        assert 0.05 < t <= 0.8

    def test_roc_returns_roc_type(self):
        _, tt = _select_threshold(
            np.array([0.1, 0.9, 0.2, 0.8]),
            np.array([0, 1, 0, 1]),
            prc_threshold=False
        )
        assert tt == "roc"

    def test_single_class_does_not_crash(self):
        scores = np.array([0.1, 0.5, 0.9])
        labels = np.zeros(3, dtype=int)
        t, tt = _select_threshold(scores, labels, prc_threshold=True)
        assert tt == "prc"
        assert isinstance(t, float)


# ---------------------------------------------------------------------------
# 8. simple criterion end-to-end
# ---------------------------------------------------------------------------

class TestSimpleCriterionEndToEnd:
    def _setup(self):
        n = 4
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = B_ref[2, 3] = 1
        W_pred = np.full((n, n), 0.05)
        W_pred[0, 1] = W_pred[1, 2] = W_pred[2, 3] = 0.8
        np.fill_diagonal(W_pred, 0.0)
        return W_pred, B_ref

    def test_recall_precision_perfect(self):
        W, B = self._setup()
        bundle = _compute_simple_metrics(W, B, 0.1, True, True)
        assert bundle["recall"]    == pytest.approx(1.0)
        assert bundle["precision"] == pytest.approx(1.0)
        assert bundle["f1"]        == pytest.approx(1.0)

    def test_jaccard_perfect(self):
        W, B = self._setup()
        assert _compute_simple_metrics(W, B, 0.1, True, True)["jaccard_similarity"] == pytest.approx(1.0)

    def test_threshold_between_noise_and_signal(self):
        W, B = self._setup()
        bundle = _compute_simple_metrics(W, B, 0.1, True, True)
        assert 0.05 < bundle["threshold"] <= 0.8

    def test_fixed_threshold_respected(self):
        W, B = self._setup()
        bundle = _compute_simple_metrics(W, B, 0.5, auto_threshold=False, prc_threshold=True)
        assert bundle["threshold"] == pytest.approx(0.5)
        assert bundle["threshold_type"] == "prc"

    def test_roc_threshold_type(self):
        W, B = self._setup()
        bundle = _compute_simple_metrics(W, B, 0.1, True, prc_threshold=False)
        assert bundle["threshold_type"] == "roc"

    def test_high_accuracy(self):
        W, B = self._setup()
        assert _compute_simple_metrics(W, B, 0.1, True, True)["accuracy"] > 0.9


# ---------------------------------------------------------------------------
# 9. all_paths criterion end-to-end
# ---------------------------------------------------------------------------

class TestAllPathsCriterionEndToEnd:
    def _chain(self):
        n = 4
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = B_ref[2, 3] = 1
        W_pred = np.full((n, n), 0.02)
        W_pred[0, 1] = W_pred[1, 2] = W_pred[2, 3] = 0.8
        np.fill_diagonal(W_pred, 0.0)
        return W_pred, B_ref

    def test_recall_one_for_perfect_chain(self):
        W, B = self._chain()
        W_reach = modified_floyd_warshall(W)
        B_reach = floyd_warshall_closure(B)
        bundle = _compute_all_paths_metrics(W, W_reach, B_reach, 0.1, True, True)
        assert bundle["recall"] == pytest.approx(1.0)

    def test_jaccard_for_perfect_chain_with_noise(self):
        """
        0.02 noise causes PRC to select threshold=0.02, predicting all entries.
        Jaccard = B_reach.sum() / n^2 = 10/16 = 0.625.
        """
        W, B = self._chain()
        W_reach = modified_floyd_warshall(W)
        B_reach = floyd_warshall_closure(B)
        bundle = _compute_all_paths_metrics(W, W_reach, B_reach, 0.1, True, True)
        n = 4
        expected_jaccard = float(B_reach.sum()) / (n * n)
        assert bundle["jaccard_similarity"] == pytest.approx(expected_jaccard, abs=1e-9)

    def test_indirect_paths_recovered(self):
        W, B = self._chain()
        B_reach = floyd_warshall_closure(B)
        assert B_reach[0, 2] == 1
        assert B_reach[0, 3] == 1


# ---------------------------------------------------------------------------
# 10. compute_graph_sim_metrics schema
# ---------------------------------------------------------------------------

class TestComputeGraphSimMetricsSchema:
    def _dfs(self):
        states = ["A", "B", "C", "D"]
        ref = pd.DataFrame(0, index=states, columns=states, dtype=int)
        ref.loc["A", "B"] = ref.loc["B", "C"] = 1
        pred = pd.DataFrame(0.05, index=states, columns=states, dtype=float)
        pred.loc["A", "B"] = 0.9; pred.loc["B", "C"] = 0.85
        return pred, ref, states

    def test_returns_single_and_multi_step(self):
        p, r, s = self._dfs()
        result = compute_graph_sim_metrics(p, r, s)
        assert "single_step" in result and "multi_step" in result

    def test_required_keys_single_step(self):
        p, r, s = self._dfs()
        required = {"threshold", "threshold_type", "accuracy", "precision",
                    "recall", "f1", "auc_roc", "auc_prc", "jaccard_similarity"}
        assert required.issubset(set(compute_graph_sim_metrics(p, r, s)["single_step"]))

    def test_required_keys_multi_step(self):
        p, r, s = self._dfs()
        required = {"threshold", "threshold_type", "accuracy", "precision",
                    "recall", "f1", "auc_roc", "auc_prc", "jaccard_similarity"}
        assert required.issubset(set(compute_graph_sim_metrics(p, r, s)["multi_step"]))

    def test_custom_criteria_simple_only(self):
        p, r, s = self._dfs()
        result = compute_graph_sim_metrics(p, r, s, threshold_criteria=["simple"])
        assert "single_step" in result and "multi_step" not in result

    def test_missing_states_zero_filled(self):
        states = ["A", "B", "C"]
        ref = pd.DataFrame(0, index=states, columns=states, dtype=int)
        ref.loc["A", "B"] = 1
        pred_partial = pd.DataFrame({"A": [0.0, 0.9], "B": [0.1, 0.0]}, index=["A", "B"])
        result = compute_graph_sim_metrics(pred_partial, ref, states)
        assert "single_step" in result and "multi_step" in result

        report = validate_predicted_label_space(pred_partial, states)
        assert report["missing_reference_nodes"] == ["C"]
        assert report["zero_filled_missing_reference_nodes"] is True

    def test_stage_labels_raise(self):
        states = ["A", "B"]
        ref = pd.DataFrame(0, index=states, columns=states, dtype=int)
        pred = pd.DataFrame(
            0.0,
            index=["stage_00_A", "stage_01_A_to_B"],
            columns=["stage_00_A", "stage_01_A_to_B"],
        )
        with pytest.raises(ValueError, match="expanded/stage labels"):
            compute_graph_sim_metrics(pred, ref, states)

    def test_prc_threshold_type_default(self):
        p, r, s = self._dfs()
        result = compute_graph_sim_metrics(p, r, s)
        assert result["single_step"]["threshold_type"] == "prc"
        assert result["multi_step"]["threshold_type"]  == "prc"


# ---------------------------------------------------------------------------
# 11. lineage_metrics.json schema from eval_lineage.run_lineage_evaluation
# ---------------------------------------------------------------------------

class TestLineageMetricsJsonSchema:
    def _run(self, tmp_path: Path) -> dict:
        graph = {
            "_meta": {"version": "test", "state_key": "state"},
            "nodes": [
                {"id": "S1", "status": "confirmed"},
                {"id": "S2", "status": "confirmed"},
                {"id": "S3", "status": "confirmed"},
            ],
            "edges": [
                {"source": "S1", "target": "S2", "weight": 1.0, "confidence": "high"},
                {"source": "S2", "target": "S3", "weight": 1.0, "confidence": "high"},
            ],
        }
        (tmp_path / "reference_graph.json").write_text(json.dumps(graph))

        states = ["S1", "S2", "S3"]
        pred = pd.DataFrame(0.05, index=states, columns=states)
        pred.loc["S1", "S2"] = 0.9; pred.loc["S2", "S3"] = 0.85
        pred.to_csv(tmp_path / "state_transition_matrix.csv")
        (tmp_path / "lineage_graph_edges.csv").write_text(
            "source_state,target_state,weight\nS1,S2,0.9\nS2,S3,0.85\n"
        )

        from benchmark.evaluation.eval_lineage import run_lineage_evaluation
        return run_lineage_evaluation(
            state_transition_matrix_path=str(tmp_path / "state_transition_matrix.csv"),
            lineage_graph_edges_path=str(tmp_path / "lineage_graph_edges.csv"),
            output_dir=str(tmp_path),
            reference_graph_path=str(tmp_path / "reference_graph.json"),
            adata=None,
            edge_confidence_mode="high_only",
            cell_state_key="state",
        )

    def test_metric_protocol(self, tmp_path):
        assert self._run(tmp_path)["metric_protocol"] == "sctimebench_graph_sim"

    def test_graph_metrics_structure(self, tmp_path):
        gm = self._run(tmp_path)["graph_metrics"]
        assert gm is not None and "single_step" in gm and "multi_step" in gm

    def test_backward_compat_aliases_present(self, tmp_path):
        m = self._run(tmp_path)
        for key in ("auroc", "auprc", "jaccard_similarity",
                    "single_step_recovery", "multi_step_recovery"):
            assert key in m

    def test_auroc_equals_single_step_auc_roc(self, tmp_path):
        m = self._run(tmp_path)
        assert m["auroc"] == pytest.approx(m["graph_metrics"]["single_step"]["auc_roc"])

    def test_auprc_equals_single_step_auc_prc(self, tmp_path):
        m = self._run(tmp_path)
        assert m["auprc"] == pytest.approx(m["graph_metrics"]["single_step"]["auc_prc"])

    def test_jaccard_equals_single_step_jaccard(self, tmp_path):
        m = self._run(tmp_path)
        assert m["jaccard_similarity"] == pytest.approx(
            m["graph_metrics"]["single_step"]["jaccard_similarity"]
        )

    def test_json_written_and_parseable(self, tmp_path):
        self._run(tmp_path)
        content = json.loads((tmp_path / "lineage_metrics.json").read_text())
        assert content["metric_protocol"] == "sctimebench_graph_sim"

    def test_legacy_topk_field_present(self, tmp_path):
        assert "legacy_jaccard_similarity_topk" in self._run(tmp_path)

    def test_no_bare_jaccard_similarity_topk(self, tmp_path):
        assert "jaccard_similarity_topk" not in self._run(tmp_path)

    def test_prediction_label_source_is_sctimebench_zero_fill(self, tmp_path):
        m = self._run(tmp_path)
        assert (
            m["prediction_label_source"]
            == "reference_graph_nodes_sctimebench_zero_filled"
        )

    def test_single_step_fields_finite(self, tmp_path):
        gm = self._run(tmp_path)["graph_metrics"]["single_step"]
        for field in ("auc_roc", "auc_prc", "jaccard_similarity", "precision", "recall", "f1"):
            val = gm.get(field)
            assert val is not None and math.isfinite(float(val)), f"single_step.{field}={val!r}"

    def test_multi_step_fields_finite(self, tmp_path):
        gm = self._run(tmp_path)["graph_metrics"]["multi_step"]
        for field in ("auc_roc", "auc_prc", "jaccard_similarity", "precision", "recall", "f1"):
            val = gm.get(field)
            assert val is not None and math.isfinite(float(val)), f"multi_step.{field}={val!r}"


# ---------------------------------------------------------------------------
# 12. align_predicted_to_reference
# ---------------------------------------------------------------------------

class TestAlignPredictedToReference:
    def test_missing_state_zero_filled(self):
        pred = pd.DataFrame({"A": [0.5, 0.3], "B": [0.2, 0.1]}, index=["A", "B"])
        W = align_predicted_to_reference(pred, ["A", "B", "C"])
        assert W.shape == (3, 3)
        assert np.allclose(W[2, :], 0.0)
        assert np.allclose(W[:, 2], 0.0)

    def test_extra_state_dropped(self):
        pred = pd.DataFrame(
            {"A": [0.9, 0.1, 0.0], "B": [0.2, 0.8, 0.0], "X": [0.0] * 3},
            index=["A", "B", "X"]
        )
        assert align_predicted_to_reference(pred, ["A", "B"]).shape == (2, 2)

    def test_ordering_matches_ref(self):
        pred = pd.DataFrame({"B": [0.7, 0.0], "A": [0.0, 0.3]}, index=["B", "A"])
        W = align_predicted_to_reference(pred, ["A", "B"])
        assert W[0, 0] == pytest.approx(0.3)
        assert W[1, 1] == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# 13. Self-loop exclusion in predicted binary graph (scTimeBench parity)
# ---------------------------------------------------------------------------

class TestSelfLoopExclusion:
    """
    scTimeBench _build_pred_graph_with_threshold() skips self-loops:
        if source_id != target_id and prob >= threshold: adjacency[src][tgt] = 1
    np.fill_diagonal(pred_binary, 0) mirrors this for both simple and all_paths.
    """

    def test_simple_diagonal_not_predicted_even_when_score_is_high(self):
        """
        W_pred diagonal = 0.99 (above any reasonable threshold).
        B_ref diagonal  = 0 (no self-loops in reference).
        After fill_diagonal, pred_binary diagonal must be 0 -> no FP from diagonal.
        """
        n = 3
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = 1

        W = np.zeros((n, n))
        W[0, 1] = W[1, 2] = 0.9
        np.fill_diagonal(W, 0.99)   # high self-transition scores

        bundle = _compute_simple_metrics(W, B_ref, 0.5, auto_threshold=False, prc_threshold=True)
        # Perfect signal on the two reference edges, diagonal excluded.
        assert bundle["jaccard_similarity"] == pytest.approx(1.0), (
            "Jaccard should be 1.0; self-loop scores must not become predicted edges."
        )
        assert bundle["precision"] == pytest.approx(1.0)
        assert bundle["recall"]    == pytest.approx(1.0)

    def test_simple_jaccard_one_with_high_diagonal(self):
        """
        Concrete 3-state case: reference A->B and B->C.
        W_pred has strong scores on those edges and high diagonal.
        At threshold 0.5 (fixed), only A->B and B->C should be predicted.
        Jaccard vs the 3-state B_ref (diagonal=0) must be 1.0.
        """
        states = ["A", "B", "C"]
        W_pred = pd.DataFrame(0.0, index=states, columns=states)
        W_pred.loc["A", "B"] = 0.9
        W_pred.loc["B", "C"] = 0.9
        W_pred.loc["A", "A"] = W_pred.loc["B", "B"] = W_pred.loc["C", "C"] = 0.95

        B_ref = pd.DataFrame(0, index=states, columns=states)
        B_ref.loc["A", "B"] = B_ref.loc["B", "C"] = 1

        result = compute_graph_sim_metrics(
            W_pred, B_ref, states,
            threshold_criteria=["simple"],
            auto_threshold=False,
            edge_threshold=0.5,
        )
        assert result["single_step"]["jaccard_similarity"] == pytest.approx(1.0), (
            "High diagonal self-transition scores must not be included as predicted edges."
        )

    def test_all_paths_direct_pred_diagonal_zero_before_closure(self):
        """
        Before floyd_warshall_closure is applied, pred_binary_direct must have
        diagonal=0 (self-loops excluded). After closure, diagonal becomes 1
        because FW adds self-reachability -- not from thresholding self-loops.
        """
        n = 3
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = 1

        W_pred = np.zeros((n, n))
        W_pred[0, 1] = W_pred[1, 2] = 0.9
        np.fill_diagonal(W_pred, 0.99)

        W_reach = modified_floyd_warshall(W_pred)
        B_reach = floyd_warshall_closure(B_ref)

        bundle = _compute_all_paths_metrics(
            W_pred, W_reach, B_reach, 0.5, auto_threshold=False, prc_threshold=True
        )

        # After closure, diagonal is 1 in pred_closed and B_reach -> always TP.
        # But self-loops were NOT included in pred_binary_direct.
        # Verify by checking recall: all B_reach positives should be covered.
        assert bundle["recall"] == pytest.approx(1.0)
        # Jaccard must not be penalized by spurious self-loop FP.
        assert bundle["jaccard_similarity"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 14. AUC diagonal zeroing (scTimeBench parity)
# ---------------------------------------------------------------------------


class TestAUCDiagonalZeroing:
    """
    scTimeBench get_threshold_roc/prc zero graph_weighted_pred_adj[i,i] = 0.0
    before calling roc_auc_score / average_precision_score.
    """

    def _chain_ref(self, n=4):
        B = np.zeros((n, n), dtype=int)
        for i in range(n - 1):
            B[i, i + 1] = 1
        return B

    def test_simple_auc_invariant_to_diagonal_value(self):
        n = 4
        B_ref = self._chain_ref(n)
        W_clean = np.zeros((n, n))
        for i in range(n - 1):
            W_clean[i, i + 1] = 0.9
        W_high_diag = W_clean.copy()
        np.fill_diagonal(W_high_diag, 0.95)
        b_clean = _compute_simple_metrics(W_clean,     B_ref, 0.1, True, True)
        b_diag  = _compute_simple_metrics(W_high_diag, B_ref, 0.1, True, True)
        assert b_clean["auc_roc"] == pytest.approx(b_diag["auc_roc"], abs=1e-9)
        assert b_clean["auc_prc"] == pytest.approx(b_diag["auc_prc"], abs=1e-9)

    def test_simple_auc_would_differ_without_diagonal_zeroing(self):
        from sklearn.metrics import roc_auc_score
        n = 4
        B_ref = self._chain_ref(n)
        W_clean = np.zeros((n, n))
        for i in range(n - 1):
            W_clean[i, i + 1] = 0.9
        W_high_diag = W_clean.copy()
        np.fill_diagonal(W_high_diag, 0.95)
        y_true = B_ref.flatten().astype(int)
        auc_clean = roc_auc_score(y_true, W_clean.flatten())
        auc_diag  = roc_auc_score(y_true, W_high_diag.flatten())
        assert auc_clean != auc_diag

    def test_all_paths_auc_invariant_to_diagonal_value(self):
        n = 4
        B_ref = self._chain_ref(n)
        W_clean = np.zeros((n, n))
        for i in range(n - 1):
            W_clean[i, i + 1] = 0.9
        W_high_diag = W_clean.copy()
        np.fill_diagonal(W_high_diag, 0.95)
        B_reach = floyd_warshall_closure(B_ref)
        b_clean = _compute_all_paths_metrics(
            W_clean, modified_floyd_warshall(W_clean), B_reach, 0.1, True, True
        )
        b_diag = _compute_all_paths_metrics(
            W_high_diag, modified_floyd_warshall(W_high_diag), B_reach, 0.1, True, True
        )
        assert b_clean["auc_roc"] == pytest.approx(b_diag["auc_roc"], abs=1e-9)
        assert b_clean["auc_prc"] == pytest.approx(b_diag["auc_prc"], abs=1e-9)


# ---------------------------------------------------------------------------
# 15. all_paths thresholds original W_pred, not W_reach
# ---------------------------------------------------------------------------

class TestAllPathsThresholdsOriginalWPred:
    """
    scTimeBench _build_pred_graph_with_threshold(weighted_adjacency, threshold)
    thresholds the original W_pred, not W_reach. W_reach is for threshold
    selection only. Mathematical note: close(W_pred >= T) == close(W_reach >= T),
    so final metrics are identical either way; this test documents the expected
    data flow and serves as a regression test.
    """

    def _setup(self):
        n = 3
        W_pred = np.zeros((n, n))
        W_pred[0, 1] = W_pred[1, 2] = 0.8
        W_reach = modified_floyd_warshall(W_pred)
        B_ref = np.zeros((n, n), dtype=int)
        B_ref[0, 1] = B_ref[1, 2] = 1
        B_reach = floyd_warshall_closure(B_ref)
        return W_pred, W_reach, B_reach

    def test_w_reach_has_indirect_not_in_w_pred(self):
        W_pred, W_reach, _ = self._setup()
        assert W_pred[0, 2] == pytest.approx(0.0)
        assert W_reach[0, 2] == pytest.approx(0.8)

    def test_result_matches_w_pred_thresholding(self):
        W_pred, W_reach, B_reach = self._setup()
        T = 0.5
        expected_direct = (W_pred >= T).astype(int)
        np.fill_diagonal(expected_direct, 0)
        expected_closed = floyd_warshall_closure(expected_direct)
        bundle = _compute_all_paths_metrics(
            W_pred, W_reach, B_reach, edge_threshold=T, auto_threshold=False, prc_threshold=True
        )
        expected_jac = _jaccard_flat(expected_closed.flatten(), B_reach.flatten())
        assert bundle["jaccard_similarity"] == pytest.approx(expected_jac, abs=1e-9)

    def test_recall_is_one_at_fixed_threshold(self):
        W_pred, W_reach, B_reach = self._setup()
        bundle = _compute_all_paths_metrics(
            W_pred, W_reach, B_reach, edge_threshold=0.5, auto_threshold=False, prc_threshold=True
        )
        assert bundle["recall"] == pytest.approx(1.0)

    def test_threshold_comes_from_w_reach_scores(self):
        W_pred, W_reach, B_reach = self._setup()
        bundle = _compute_all_paths_metrics(
            W_pred, W_reach, B_reach, edge_threshold=0.1, auto_threshold=True, prc_threshold=True
        )
        thresh = bundle["threshold"]
        unique_w_reach = np.unique(W_reach.flatten())
        assert any(abs(thresh - v) < 1e-9 for v in unique_w_reach), (
            f"Threshold {thresh} should come from W_reach scores."
        )


if __name__ == "__main__":
    import sys
    import pytest as _pytest
    sys.exit(_pytest.main([__file__, "-v"]))
