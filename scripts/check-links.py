#!/usr/bin/env python3

import aiohttp
import asyncio
from pathlib import Path
import logging
from typing import List, Tuple
import time
from scripts.progress_manager import ProgressManager, TaskProgress

class LinkValidator:
    """Validates external resource links with proper error handling"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.progress = ProgressManager(logger)
        
    async def _validate_link(self, session: aiohttp.ClientSession, url: str) -> Tuple[str, int]:
        """Validate a single URL with retry logic"""
        max_retries = 3
        retry_delay = 1
        
        for attempt in range(max_retries):
            try:
                async with session.head(url, allow_redirects=True, timeout=30) as response:
                    if response.status == 200:
                        return url, 0
                    self.logger.warning(f"URL {url} returned status {response.status}")
                    
                # Retry with GET if HEAD fails
                async with session.get(url, timeout=30) as response:
                    return url, 0 if response.status == 200 else 1
                    
            except asyncio.TimeoutError:
                self.logger.warning(f"Timeout accessing {url}, attempt {attempt + 1}/{max_retries}")
                if attempt < max_retries - 1:
                    await asyncio.sleep(retry_delay)
                    retry_delay *= 2  # Exponential backoff
                    continue
                return url, 1
                
            except Exception as e:
                self.logger.error(f"Error validating {url}: {str(e)}")
                return url, 1
                
        return url, 1

    async def validate_links(self, links: List[str]) -> List[Tuple[str, int]]:
        """Validate multiple links concurrently"""
        task = TaskProgress(name="Validating external resources", total=len(links))
        self.progress.start_task(task)
        
        async with aiohttp.ClientSession() as session:
            tasks = [self._validate_link(session, url) for url in links]
            results = []
            
            for future in asyncio.as_completed(tasks):
                result = await future
                results.append(result)
                self.progress.update()
                
            self.progress.complete_task()
            return results

def main():
    """Main entry point"""
    logger = logging.getLogger("flync.links")
    validator = LinkValidator(logger)
    
    # Load required links
    appdir = Path(__file__).resolve().parent.parent
    req_links = Path(appdir) / "static" / "required_links.txt"
    
    try:
        with open(req_links, "r") as f:
            links = [line.strip() for line in f if line.strip()]
            
        if not links:
            logger.error("No links found in required_links.txt")
            raise ValueError("Empty links file")
            
        # Run validation
        results = asyncio.run(validator.validate_links(links))
        
        # Process results
        failed = [(url, code) for url, code in results if code == 1]
        if failed:
            for url, _ in failed:
                logger.error(f"URL not reachable: {url}")
            raise RuntimeError("One or more required links are not reachable")
            
        logger.info("All required links are reachable")
        
    except Exception as e:
        logger.error(f"Link validation failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()