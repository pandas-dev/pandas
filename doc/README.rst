.. _contributing.docs:

문서작업에 기여하세요!
=================================

비록 당신이 개발자가 아니어도, 문서작업에 기여하는 것만으로도 굉장히 가치있는 일입니다.
*pandas* 의 전문가가 되지 않더라도 할 수 있는 일입니다!
문서를 명확하게 만들기 위해 어떤 짧은 문구만 고침으로써 간단하고 효과적으로 기여할 수 있습니다.
다음 사람이 그 구절을 읽은만큼 당신의 기여도는 빛을 발하게 됩니다!

실제로 전문가들이 썼음에도 읽기 쉽지않은 문서들이 있습니다.
만약 당신이 이해가 안 가는 문서가 있으면,
그 부분을 당신이 확인한 후 업데이트하는 것은 다음 사람에게 도움이 될 확실한 방법입니다.

.. contents:: Table of contents:
   :local:


Pandas 문서 정보
------------------------------

이 문서는 **reStructuredText** 로 작성되어 있습니다. 이것은 보통 영어로 작성하고 `Sphinx <http://sphinx.pocoo.org/>`_ 를 사용하여 작성한 것과 같습니다. Sphinx Documentation에는 reST에 대한 훌륭한 `소개 <http://sphinx.pocoo.org/rest.html>`_ 가 있습니다. Sphinx 문서를 검토하여 문서에 대한 보다 복잡한 변경을 수행할 수 있습니다.


문서에 대해 알아야 할 다른 중요한 사항 :

- 판다 문서는 코드 자체의 문서화 문자열과 이 폴더의 문서의 두 부분으로 구성됩니다 ``pandas/doc/``.
  docstring은 개별 기능의 사용법에 대한 명확한 설명을 제공하지만, 이 폴더의 설명서는 주제별 자습서와 유사한 개요(새로운 기능, 설치 등)와 함께 제공됩니다.

- 문서화 문자열은 과학범용 파이썬 커뮤니티에서 널리 사용되는 **Numpy Docstring 표준** 을 따르고 있습니다. 이 표준은 문서화 문자열의 여러 섹션의 형식을 지정합니다. 자세한 `설명 <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ 을 보려면 이 문서 를 참조 하거나 비슷한 방식으로 기존 기능을 확장하십시오.

- 이 튜토리얼 에서는 `ipython 지시문 <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx 확장을 많이 사용합니다.   이 지시문을 사용하면 실행되는 문서에 코드를 넣을 수 있습니다.

  문서 작성 중에 예를 들면 ::
  
      .. ipython:: python
          x = 2
          x**3
          
  은 다음과 같이 랜더링 됩니다.::
  
      In [1]: x = 2
      In [2]: x**3
      Out[2]: 8

  즉, 문서 작성 중에 거의 모든 코드 예제가 실행되고 출력이 저장됩니다. 이렇게 하면 항상 최신 상태가 유지되지만 문서 작성은 좀 더 복잡해집니다.

Pandas 문서 작성 방법
-------------------------------------

요구 사항
^^^^^^^^^^^^^^

팬더의 문서를 만들려면 몇 가지 추가 요구 사항이 있습니다 : ``sphinx`` 와 ``ipython`` 을 설치해야 합니다.
`numpydoc <https://github.com/numpy/numpydoc>`_ 은 Numpy Docstring Standard (위 참조)를 따르는 문서화 문자열을 구문 분석하는 데 사용되지만 “numpydoc“ 판다 소스 코드에 로컬 사본 이 포함되어 있기 때문에 설치하지 않아도 됩니다. `nbsphinx <https://nbsphinx.readthedocs.io/>`_ 는 Jupyter 노트북을 변환하는 데 사용됩니다. 이 문서에 포함 된 노트북을 수정하려면 이 프로그램을 설치해야합니다.

또한, `선택적 종속성 <http://pandas.pydata.org/pandas-docs/dev/install.html#optional-dependencies>`_ 을 모두 설치하는 것이 좋습니다. 이이것은 필요하지 않지만 몇 가지 오류 메시지가 나타납니다. 문서의 모든 코드가 문서 빌드 중에 실행되므로 이 선택적 종속성을 사용하는 예제는 오류를 생성합니다. ``pd.show_versions ()`` 을 실행하면 모든 종속성의 설치된 버전에 대한 개요를 볼 수 있습니다.

.. 경고::

   Sphinx 버전> = 1.2.2 또는 그 이전 1.1.3이 필요합니다.

pandas 만들기
^^^^^^^^^^^^^^^^^^

환경 설정 방법에 대한 단계별 개요, pandas 코드 및 git에 대한 작업은 `개발자 페이지 <http://pandas.pydata.org/developers.html#working-with-the-code>`_ 를 참조하십시오. 
일부 문서 작업을 시작할 때는 코드를 최신 개발 버전('master')으로 업데이트해야 합니다 ::

    git fetch upstream
    git rebase upstream / master

업데이트 후에 C 확장을 다시 빌드해야하는 경우가 종종 있습니다 ::

    python setup.py build_ext —inplace

문서 작성
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

그렇다면 어떻게 문서를 만들 수 있습니까? 
콘솔에서 ``pandas / doc /`` 로컬 폴더 디렉토리로 이동하여 다음을 실행하십시오::

    python make.py html

그리고 나서 폴더에서 html 출력을 찾을 수 있습니다 ``pandas / doc / build / html /``.

그것이 문서의 모든 코드 예제를 실행하고 생성된 모든 docstring 페이지를 빌드해야하기 때문에 처음에는 꽤 오래 걸릴 것입니다. 후속 연상에서, 스핑크스는 수정 된 페이지 만 만들려고 시도합니다.

완전히 깔끔하게 빌드하려면 다음을 수행하십시오. ::

    python make.py clean
    python make.py 빌드
    
    
0.13.1부터는 ``make.py`` 문서 섹션을 하나만 컴파일하면 변경 사항을 확인하는 데 소요되는 시간을 크게 줄일 수 있습니다. 마지막 커밋 된 버전은 항상 git에서 복원 할 수 있기 때문에 불필요한 .rst 파일을 삭제하라는 메시지가 표시됩니다.

::

    #omit autosummary 및 API 섹션
    python make.py clean
    python make.py --no-api
    
    # 문서를 하나만 컴파일
    # section, indexing.rst에있는 섹션
    python make.py clean
    python make.py - 단일 지수 연동
    
비교를 위해 전체 문서 빌드에 10 분이 소요될 수 있습니다. ``-no-api`` 빌드 3 분 정도 걸릴 수 있으며, 단일 섹션 15 초가 걸릴 수 있습니다.

어디서부터 시작해야 할까?
---------------

처음 시작할 수 있는 `Docs <https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open>`_ 와 `최초의 PR로 좋은 <https://github.com/pandas-dev/pandas/issues?labels=Good+as+first+PR&sort=updated&state=open>`_ 여러 가지 문제가 있는 리스트가 있습니다.

아니면 판다를 사용하고, 문서에서 뭔가를 찾고 '이것이 향상 될 수 있습니다'라고 생각하면 그 생각에 대한 행동을 하십시오!

`메일 링리스트 <https://groups.google.com/forum/?fromgroups#!forum/pydata>`_ 에 질문 하거나 Github에 문제를 제출하십시오.
