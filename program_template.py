from __future__ import division, print_function
'''
Standard Program Template for Amber Programs

Similar to the one in CCTBX, but with some cosmetic changes
'''

import libtbx.program_template

from libtbx import citations

# =============================================================================
class ProgramTemplate(libtbx.program_template.ProgramTemplate):

  epilog = '''
For additional help, you can contact the developers at help@phenix-online.org
'''

  # ---------------------------------------------------------------------------
  @staticmethod
  def show_template_citation(text_width=80,
                             logger=None,
                             citation_format='default'):
    assert(logger is not None)

    print('\nGeneral citation for Amber:', file=logger)
    print('-'*text_width, file=logger)
    print('', file=logger)
    citations.show_citation(citations.citations_db['amber'],
                            out=logger,
                            format=citation_format)
    libtbx.program_template.ProgramTemplate.show_template_citation(
      text_width=text_width, logger=logger, citation_format=citation_format)

# =============================================================================
