"""Synthetic inputs shared by the offline unit suite."""

import importlib.util
import sys
from pathlib import Path

import pytest
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks" / "workflow"))
sys.path.insert(0, str(ROOT / "notebooks"))
import workflow_utils as workflow_utils
import cs3_utils as cs3_utils

spec = importlib.util.spec_from_file_location(
    "prepare_events", ROOT / "downloads" / "prepare_events.py"
)
prepare_events = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_events)


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    def reject(*args, **kwargs):
        raise AssertionError("Unit tests must not contact a remote service")

    monkeypatch.setattr(requests.sessions.Session, "request", reject)


@pytest.fixture
def events():
    def make(**changes):
        row = dict(
            ID="first",
            GeneID='"ENSG00000000001.2"',
            geneSymbol="GENE",
            chr="chr1",
            strand="+",
            exonStart_0base=100,
            exonEnd=120,
            upstreamES=10,
            upstreamEE=30,
            downstreamES=200,
            downstreamEE=220,
            IJC_SAMPLE_1="10,20",
            SJC_SAMPLE_1="2,3",
            IJC_SAMPLE_2="5,5",
            SJC_SAMPLE_2="10,10",
            IncLevel1="0.7,0.9",
            IncLevel2="0.4,0.6",
            IncLevelDifference=0.2,
            FDR=0.01,
            PValue=0.01,
        )
        row.update(changes)
        return row

    return make


class OffsetLift:
    def convert_coordinate(self, chrom, pos):
        return [(chrom, pos + 1000, "+", 100)]


@pytest.fixture
def lift():
    return OffsetLift()


@pytest.fixture
def locus():
    exons = {k: list(v["grch38"][2:]) for k, v in workflow_utils.EXONS.items()}
    a, b = exons["exon13"], exons["exon24nt"]
    outside = [[85970000, 85970100], [85996000, 85996100]]
    tx = {
        "include": dict(
            name="PICALM-1",
            type="protein_coding",
            tags=["MANE_Select"],
            exons=sorted(outside + [a, b]),
        ),
        "skip": dict(name="PICALM-2", type="protein_coding", tags=[], exons=outside),
        "short": dict(
            name="PICALM-3",
            type="protein_coding",
            tags=[],
            exons=sorted(outside + [[a[0], a[1] - 21]]),
        ),
        "noncoding": dict(
            name="PICALM-4", type="retained_intron", tags=[], exons=outside
        ),
    }
    return dict(
        exons=exons,
        transcripts=tx,
        windows=[
            ["QKI_K562_5007", a[1], a[1] + 100, 4, 0.01],
            ["QKI_HepG2_ENCSR570WLM", b[0] - 100, b[0], 2.5, 0.01],
        ],
    )


@pytest.fixture
def knockdown():
    st = dict(
        psi_ctl=[0.3, 0.4],
        psi_kd=[0.6, 0.7],
        IJC_ctl=[20, 30],
        SJC_ctl=[50, 60],
        IJC_kd=[60, 70],
        SJC_kd=[20, 30],
        dpsi=0.3,
        fdr=0.01,
        pvalue=0.001,
    )
    return {
        k: dict(log2fc=-2.0, exon13=dict(st), exon24nt=dict(st, dpsi=-0.3, fdr=0.6))
        for k in workflow_utils.KNOCKDOWNS
    }
