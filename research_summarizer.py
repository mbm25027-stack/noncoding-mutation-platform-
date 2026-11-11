"""
Production-Grade Research Paper Summarizer & Knowledge Base
----------------------------------------------------------
1. Ingestion: PDF parsing with section splitting (PyMuPDF / fitz)
2. Schema Extraction: Forced Pydantic JSON structure using OpenAI API
3. Verification: Metric cross-checking against source document
4. Vector Storage: Multi-paper search & indexing foundation
"""

import os
import re
from typing import List, Dict, Optional
from pydantic import BaseModel, Field
import fitz  # PyMuPDF
from openai import OpenAI


# 1. STRUCTURED SCHEMA DEFINITION
class QuantitativeMetric(BaseModel):
    metric_name: str = Field(description="Name of benchmark or metric (e.g., Accuracy, F1 Score, Latency)")
    value: str = Field(description="Exact reported value or result with units")
    context: str = Field(description="Contextual comparison (e.g., Outperformed baseline by 3.2%)")

class PaperSummarySchema(BaseModel):
    title: str = Field(description="Title of the research paper")
    authors: List[str] = Field(description="List of primary authors")
    core_problem: str = Field(description="2-3 sentence summary of the core research problem")
    key_contributions: List[str] = Field(description="List of novel technical contributions")
    methodology_overview: str = Field(description="Concise description of the approach or architecture")
    quantitative_results: List[QuantitativeMetric] = Field(description="Key quantitative benchmarks and results")
    limitations_and_drawbacks: List[str] = Field(description="Stated limitations, assumptions, or failure modes")
    practical_takeaways: List[str] = Field(description="Key takeaways for practitioners and future work")


# 2. PDF INGESTION & SECTION PARSING
class PDFIngester:
    @staticmethod
    def extract_sections(pdf_path: str) -> Dict[str, str]:
        """Extracts text from a PDF and splits it into logical sections."""
        doc = fitz.open(pdf_path)
        full_text = [page.get_text() for page in doc]
        combined_text = "\n".join(full_text)
        
        section_headings = ["ABSTRACT", "INTRODUCTION", "METHOD", "METHODOLOGY", "EXPERIMENTAL SETUP", "RESULTS", "DISCUSSION", "CONCLUSION"]
        pattern = r'(?i)\n(?=(' + '|'.join(section_headings) + r')\b)'
        splits = re.split(pattern, combined_text)
        
        sections = {}
        current_title = "HEADER/ABSTRACT"
        for chunk in splits:
            clean_chunk = chunk.strip()
            if not clean_chunk:
                continue
            first_line = clean_chunk.split('\n')[0].upper()
            matched_heading = next((h for h in section_headings if h in first_line), None)
            if matched_heading:
                current_title = matched_heading
            sections[current_title] = sections.get(current_title, "") + "\n" + clean_chunk
                
        return sections if sections else {"FULL_DOCUMENT": combined_text}


# 3. STRUCTURED LLM SUMMARIZER
class ResearchSummarizerPipeline:
    def __init__(self, api_key: Optional[str] = None):
        self.client = OpenAI(api_key=api_key or os.environ.get("OPENAI_API_KEY"))

    def extract_structured_summary(self, raw_text: str) -> PaperSummarySchema:
        """Extracts strictly-typed structured paper insights using Pydantic."""
        system_prompt = (
            "You are an expert senior researcher. Extract precise, rigorous technical details "
            "from the provided paper text. Do not invent details; extract only facts grounded in the text."
        )
        
        completion = self.client.beta.chat.completions.parse(
            model="gpt-4o-2024-08-06",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": f"Extract structured insights from this paper text:\n\n{raw_text[:25000]}"}
            ],
            response_format=PaperSummarySchema,
        )
        return completion.choices[0].message.parsed

    # 4. FACT & METRIC VERIFICATION PASS
    def verify_metrics(self, summary: PaperSummarySchema, original_text: str) -> List[Dict[str, str]]:
        """Cross-verifies extracted numerical metrics back against original source text."""
        verification_results = []
        for metric in summary.quantitative_results:
            found = metric.value.lower() in original_text.lower() or any(num in original_text for num in re.findall(r'\d+\.?\d*', metric.value))
            verification_results.append({
                "metric": metric.metric_name,
                "value": metric.value,
                "verified": found,
                "status": "Verified in source text" if found else "Flagged for manual review"
            })
        return verification_results


# 5. KNOWLEDGE BASE & SEARCH INDEX
class PaperKnowledgeBase:
    def __init__(self):
        self.papers: List[PaperSummarySchema] = []

    def add_paper(self, paper_summary: PaperSummarySchema):
        self.papers.append(paper_summary)

    def search(self, query: str) -> List[PaperSummarySchema]:
        """Finds indexed papers matching key terms across core problem and methodology."""
        q = query.lower()
        return [
            p for p in self.papers 
            if q in p.core_problem.lower() 
            or q in p.methodology_overview.lower()
            or any(q in c.lower() for c in p.key_contributions)
        ]


# USAGE EXAMPLE
if __name__ == "__main__":
    # Set your API key: export OPENAI_API_KEY="your-key"
    pipeline = ResearchSummarizerPipeline()
    kb = PaperKnowledgeBase()

    # Process a PDF paper
    # pdf_path = "path/to/paper.pdf"
    # sections = PDFIngester.extract_sections(pdf_path)
    # full_text = "\n".join(sections.values())
    # summary = pipeline.extract_structured_summary(full_text)
    # verification = pipeline.verify_metrics(summary, full_text)
    # kb.add_paper(summary)

    print("Pipeline ready. Execute pipeline on any PDF paper path.")