.. _contributing.docs:

문서작업에 기여하세요!
=================================

만약 당신이 개발자가 아니라면, 문서작업에 기여하는 것만으로도 굉장히 가치있는 일입니다.
“pandas”의 전문가가 되지 않더라도 할 수 있는 일입니다!
문서를 명확하게 만들기 위해 어떤 짧은 문구만 고침으로써 간단하고 효과적으로 기여할 수 있습니다.
다음 사람이 그 구절을 읽은만큼 당신의 기여도는 빛을 발하게 됩니다!

실제로 전문가들이 썼음에도 읽기 쉽지않은 문서들이 있습니다.
만약 당신이 이해가 안 가는 문서가 있으면,
그 부분을 당신이 확인한 후 업데이트하는 것은 다음 사람에게 도움이 될 확실한 방법입니다.

.. contents:: Table of contents:
   :local:


팬더 문서를 보려면
------------------------------

문서는 ** reStructuredText **로 작성되어 있습니다.
보통 영어로 쓰여진`Sphinx <http://sphinx.pocoo.org/>`__을 사용해 구축되어 있습니다. 더
Sphinx Documentation는 뛰어난`reST 소개
<http://sphinx.pocoo.org/rest.html>`__. Sphinx의 문서를 검토하고 더 많은 기능을 수행하는
복잡한 문서의 변경도 마찬가지입니다.

문서에 대해 알아야 할 다른 중요한 점 :

 - 팬더 문서는 두 부분으로 구성되어 있습니다. 코드의 문서 문자열
  이 폴더의``pandas / doc /``안의 문서와 그 문서를 삭제합니다.

  문서 문자열은 개인의 사용 방법을 명확하게 설명하고 있습니다
  이 폴더의 문서는 튜토리얼 같다
  주제별 개요 및 기타 정보 (새로운 기능
  설치 등).

- 문서화 문자열은 과학범용 파이썬 커뮤니티에서 널리쓰이는 **Numpy Docstring Standard**를 따르고 있습니다.
  이 표준은 각기 다른 부분의 문서화 문자열 양식을 구체화하고 있습니다.
  이 문서를 보시면
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  자세한 설명이나 간단히 적용하고 싶은 함수들을 볼 수 있습니다.

- 이 튜토리얼에서는`ipython 지시문
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx 확장자.
  이 지시문을 사용하면 실행되는 문서에 코드를 넣을 수 있습니다
  문서 작성 중에 예를 들면 :

  ::

      .. ipython:: python

          x = 2
          x**3

  will be rendered as

  ::

      In [1]: x = 2

      In [2]: x**3
      Out[2]: 8

이 문서의 거의 모든 코드 예제가 항상 실행되는 것을 의미합니다.
  출력이 저장됩니다). 이렇게하여 그들은 항상 최신 상태로됩니다.
  하지만 독 빌딩을 좀 더 복잡합니다.


팬더의 문서를 작성하는 방법
-------------------------------------

요구 사항
^^^^^^^^^^^^^^

팬더의 문서를 만들려면 몇 가지 추가 요구 사항이 있습니다.
``sphinx``와``ipython``가 설치되어 있습니다. `numpydoc
<https://github.com/numpy/numpydoc>`_ 그 문서 문자열을 해석하는 데 사용됩니다
Numpy Docstring Standard (위 참조)에 따르십시오. 그러나 설치할 필요가 없습니다
이것은 "numpydoc"의 로컬 복사본이 팬더 소스에 포함되어 있기 때문입니다
코드. `nbsphinx <https://nbsphinx.readthedocs.io/>`_ 변환에 사용됩니다
목성 노트. 중 하나를 변경하려면 설치해야합니다.
그 노트에 포함되어 있습니다.

또한, 모든 옵션의 종속성을 가지는 것을 권장합니다
<http://pandas.pydata.org/pandas-docs/dev/install.html#optional-dependencies>`_

설치되어 있습니다. 이것은 필요는 없지만, 어떤 오류가 표시되는 것에주의하십시오
메시지. 문서의 모든 코드가 doc에서 실행되기 때문에
이 옵션의 종속성을 사용하는 예는 오류를 생성합니다.
``pd.show_versions ()``을 실행하면 설치되어있는 모든 버전의 개요가 표시됩니다
의존성.

.. 경고 ::

   Sphinx 버전> = 1.2.2 또는 그 이전 1.1.3이 필요합니다.

pandas building
^^^^^^^^^^^^^^^^^^

환경을 설정하는 방법, 작업 방법에 대한 단계별 개요는
팬더의 코드와 git는`개발자 페이지
<http://pandas.pydata.org/developers.html#working-with-the-code>`_.
일부 문서에서 작업을 시작할 때 코드를 최신 버전으로 업데이트하십시오
개발 버전 ( 'master') ::

    git fetch upstream
    git rebase upstream / master

자주 업데이트 후 C 확장을 재구성해야합니다 :

    python setup.py build_ext --inplace

문서 작성
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

그럼 어떻게 문서를 작성하고 있습니까? 당신의 로컬 폴더로 이동
콘솔의``pandas / doc /``디렉토리로 이동합니다 ::

    python make.py html

그리고 html 출력은``pandas / doc / build / html /``폴더에 있습니다.

모든 코드를 실행해야하므로 처음에는 꽤 시간이 걸립니다
문서의 예를 참조하여 생성 된 모든 docstring 페이지를 구축합니다.
다음 질문은 스핑크스
변경되었습니다.

당신이 완전히 깨끗한 빌드하고 싶다면 ::

    python make.py clean
    python make.py 빌드

0.13.1에서``make.py``에 단일 섹션 만 컴파일하도록 지시 할 수 있습니다
변경 내용을 확인하기위한 처리 시간을 크게 단축합니다.
필요하지 않은`.rst` 파일을 삭제하도록 요청합니다.
마지막으로 커밋 된 버전은 항상 git에서 복원 할 수 있습니다.

::

    #omit autosummary 및 API 섹션
    python make.py clean
    python make.py --no-api

    # 문서를 하나만 컴파일
    # section, indexing.rst에있는 섹션
    python make.py clean
    python make.py - 단일 지수 연동

비교를 위해 전체 문서 빌드 10 분 정도 걸립니다. ``-no-api`` 빌드
3 분 정도 걸릴 수 있으며, 하나의 섹션에 15 초 정도 걸립니다.

어디서 시작해야할까?
---------------


`Docs에는 몇 가지 문제가 있습니다
<https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open>`_
최초의 PR로 좋은
<https://github.com/pandas-dev/pandas/issues?labels=Good+as+first+PR&sort=updated&state=open>`_
당신은 어디에서 시작할 수 있습니까?

또는 당신은 자신의 생각을 가지고 있을지도 모릅니다. 뭔가를 찾고 팬더를 사용하여
문서에서 "이것은 개선 할 수있다"라고 생각하고, 어떤 일을하자
그것에 대해!

`메일 링리스트에 대한 질문
<https://groups.google.com/forum/?fromgroups#!forum/pydata>`_ 또는
Github 문제
